#pragma once
#include <concepts>
#include <cassert>
#include <cstdint>
#include <cstddef>
#include <chrono>

#define func constexpr auto [[nodiscard]]
#define let const auto
#define mut auto
#define cast(...) (__VA_ARGS__)


func min(std::uint64_t left, std::uint64_t right)
{
    return left < right ? left : right;
}

func max(std::uint64_t left, std::uint64_t right)
{
    return left > right ? left : right;
}

func min(double left, double right)
{
    return left < right ? left : right;
}

func max(double left, double right)
{
    return left > right ? left : right;
}

template <typename Fn, typename To = void>
concept thunk = requires(Fn fn)
{
    { fn() } -> std::convertible_to<To>;
};

namespace detail
{
    using u64 = std::uint64_t;
    using f64 = double;
    using time_point = std::chrono::steady_clock::time_point;
    using ms = std::chrono::milliseconds;
    using ns = std::chrono::nanoseconds;

    func clock() -> time_point
    {
        if(!std::is_constant_evaluated())
            return std::chrono::high_resolution_clock::now();
        else 
            return time_point{};
    }

    func ns_diff(time_point from, time_point to) -> std::uint64_t
    { 
        if(!std::is_constant_evaluated())
            return duration_cast<ns>(to - from).count(); 
        else
            return 0;
    };

    enum class Benchmark_Stage
    {
        Initial,
        Blocked,
    };
}

template <thunk Fn>
func ellapsed_time(Fn fn) -> detail::ns
{
    let from = detail::clock();
    fn();
    return detail::clock() - from;
}

template <thunk Fn>
func count_iters(detail::ms ms, Fn fn, detail::ms warm_up_ms = ms / 10, detail::ms base_block_size = 1, detail::u64 optimal_samples = 5) -> std::uint64_t
{
    using namespace detail;
    using detail::Benchmark_Stage;

    ms *= 1000'000;
    warm_up_ms *= 1000'000;

    u64 total_ms = ms + warm_up_ms;
    u64 iters = 0;
    let start = detail::clock();
    mut to_time = warm_up_ms;
    u64 block_size = base_block_size;

    Benchmark_Stage stage = Benchmark_Stage::Initial;
    assert(optimal_samples > 0);

    if(warm_up_ms == 0)
    {
        stage = Benchmark_Stage::Blocked;
        to_time = total_ms;
    }

    while(true)
    {
        for(u64 i = 0; i < block_size; i++)
            fn();

        let ellapsed = detail::ns_diff(start, clock());
        iters += block_size;
        if(ellapsed > to_time)
        {
            if(stage == Benchmark_Stage::Initial)
            {
                block_size = iters * (total_ms - ellapsed) / (ellapsed * optimal_samples);
                block_size = max(block_size, 1);
                iters = 0; //discard so far
                to_time = total_ms;
                stage = Benchmark_Stage::Blocked;
            }
            else
                break;
        }
    }

    return iters;
}

struct Benchmark_Constants
{
    std::uint64_t optimal_samples = 100;
    std::uint64_t base_block_size = 10;
    double base_history = 0.9;
    double base_delta = 0.01;
    double history_softening_grow = 0.3;
    double history_hardening_grow = 0.3;
    double delta_softening_grow = 1.05;
    double delta_hardening_grow = 0.90;
    double history_min = 0.001;
    double history_max = 0.999;
    double delta_max = 0.15;
};

struct Benchmark_Stats
{
    double history_size;
    double avg_history_size;
    double max_history_size;
    double min_history_size;

    double delta;
    double avg_delta;
    double max_delta;
    double min_delta;
};

struct Benchmark_Results
{
    std::uint64_t total_real_time;  
    std::uint64_t total_time;  
    std::uint64_t total_iters;  

    std::uint64_t block_size;
    std::uint64_t total_samples;
    std::uint64_t accepted_samples;

    std::uint64_t iters;
    std::uint64_t time;

    double iters_per_ms;
    double ms_per_iter;

    double iters_per_ns;
    double ns_per_iter;
};

struct Benchmark_Combined_Results
{
    Benchmark_Results results;
    Benchmark_Stats stats;
};


template <bool DO_STATS = false, bool DO_HISTORY_CHANGE = true, bool DO_HISTORY_SIZE = true, thunk Fn = int(*)()>
func benchmark(std::uint64_t optimal_time, std::uint64_t min_time, std::uint64_t max_time, std::uint64_t optimal_samples, Fn fn, Benchmark_Constants constants) -> std::conditional_t<DO_STATS, Benchmark_Combined_Results, Benchmark_Results>
{
    using u64 = std::uint64_t;
    using f64 = double;
    using namespace detail;
    using detail::Benchmark_Stage;

    min_time *= 1000'000;
    max_time *= 1000'000;
    optimal_time *= 1000'000;

    let sum_geometrict_series = [](f64 damping){
        return 1.0 / (1.0 - damping);
    };

    //constexpr bool DO_HISTORY_CHANGE = true;
    constexpr bool DO_DELTA_CHANGE = true;

    f64 required_delta = constants.base_delta;
    f64 history_damping = constants.base_history;
    f64 history_size = sum_geometrict_series(history_damping);

    f64 history_size_highwater_mark = 0.0;
    f64 history_size_lowwater_mark = 999999;
    f64 history_size_running_sum = 0;
    f64 delta_highwater_mark = 0.0;
    f64 delta_lowwater_mark = 999999;
    f64 delta_running_sum = 0;

    u64 block_size = 1;
    u64 total_iters = 0;
    u64 total_time = 0;
    u64 accepted_samples = 0;
    u64 total_samples = 0;
    u64 last_iters = 0;
    u64 to_time = min_time;
    Benchmark_Stage stage = Benchmark_Stage::Initial;

    u64 reported_iters = 0;
    u64 reported_time = 0;
    f64 running_average_delta = 0;
    f64 prev_time = 0.0;

    let start = detail::clock();
    mut last = start;

    //Algorhitm:
    // (1) Run the given function untill min_time: This will lets us calculate block size 
    //     (ammount of times to run the function in sequence counting as a single sample)
    // 
    // (2) Run block and check its runtime difference against previous run. 
    //     If the difference is under some delta we add the sample and make the delta smaller
    //     If the difference is over delta we ignore the sample and increase delta 
    //     We also add the curren difference into the running_average_delta and damp it by history_damping
    //       so that we have sample to check the next runs against.
    //     We also vary ammount of damping of the running_average_delta (history) which lets us control how
    //       far into the past we are considering results. This is important for functions with big variations in runtime.
    //     
    // (3) We report back the obtained data and statistics

    while(true)
    {
        for(u64 i = 0; i < block_size; i++)
            fn();

        total_iters += block_size;

        let now = detail::clock();
        const u64 ellapsed = detail::ns_diff(start, now);
        if(ellapsed > to_time)
        {
            if(ellapsed > max_time)
                break;

            if(stage == Benchmark_Stage::Initial)
            {
                total_time += ellapsed;
                if(optimal_time > ellapsed)
                {
                    block_size = total_iters * (optimal_time - ellapsed) / (ellapsed * optimal_samples);
                    block_size = max(block_size, 1);
                }
                else
                    block_size = 1;

                to_time = 0;
                stage = Benchmark_Stage::Blocked;
            }
            else if(stage == Benchmark_Stage::Blocked)
            {   
                //calculate the current block stats
                const u64 current_iter_time = detail::ns_diff(last, now);
                const u64 current_iter_count = total_iters - last_iters;

                //calculate delta squared against the last block
                const f64 curr_time = cast(f64) current_iter_time / current_iter_count;
                const f64 delta = curr_time - prev_time;
                const f64 delta2 = delta * delta;
                const f64 curr_time2 = curr_time * curr_time;

                running_average_delta *= history_damping;
                running_average_delta += delta2;


                //if delta is in required range the block is added else is ignored
                f64 relative_average_delta;

                if constexpr(DO_HISTORY_SIZE)
                    relative_average_delta = running_average_delta / (curr_time2 * history_size);
                else
                    relative_average_delta = running_average_delta / curr_time2;

                if(relative_average_delta < required_delta)
                {
                    reported_iters += current_iter_count;
                    reported_time += current_iter_time;
                    accepted_samples++;

                    //harden requirements for delta and increase history size => less blocks will get accepted
                    if constexpr (DO_DELTA_CHANGE)
                        required_delta *= constants.delta_hardening_grow;;
                    if constexpr (DO_HISTORY_CHANGE)
                        history_damping += (constants.history_min - history_damping) * constants.history_hardening_grow;
                }
                else
                {
                    //soften requirements for delta and decrease history size => more blocks will get accepted
                    if constexpr (DO_DELTA_CHANGE)
                        required_delta *= constants.delta_softening_grow;;
                    if constexpr (DO_HISTORY_CHANGE)
                        history_damping += (constants.history_max - history_damping) * constants.history_softening_grow;
                }

                //update constants and track stats
                if constexpr (DO_HISTORY_CHANGE)
                    history_size = sum_geometrict_series(history_damping);

                if constexpr (DO_HISTORY_CHANGE && DO_STATS)
                {
                    history_size_running_sum += history_size;
                    history_size_highwater_mark = max(history_size_highwater_mark, history_size);
                    history_size_lowwater_mark = min(history_size_lowwater_mark, history_size);
                }

                if constexpr (DO_DELTA_CHANGE && DO_STATS)
                {
                    delta_running_sum += required_delta;
                    delta_highwater_mark = max(delta_highwater_mark, required_delta);
                    delta_lowwater_mark = min(delta_lowwater_mark, required_delta);
                }

                prev_time = curr_time;
                total_samples++;
                total_time += current_iter_time;
            }

            last = detail::clock();
            last_iters = total_iters;
        }
    }

    if(reported_iters == 0)
    {
        reported_iters = total_iters;
        reported_time = detail::ns_diff(start, detail::clock());
    }

    const f64 ns_per_iter = cast(f64) reported_iters / cast(f64) reported_time;
    const f64 ms_per_iter = ns_per_iter / 1'000'000;
    const u64 total_real_time = detail::ns_diff(start, detail::clock());

    let reuslts = Benchmark_Results{

        .total_real_time = total_real_time,  
        .total_time = total_time,  
        .total_iters = total_iters,  

        .block_size = block_size,
        .total_samples = total_samples,
        .accepted_samples = accepted_samples,

        .iters = reported_iters,
        .time = reported_time,

        .iters_per_ms = 1 / ms_per_iter,
        .ms_per_iter = ms_per_iter,

        .iters_per_ns = 1 / ns_per_iter,
        .ns_per_iter = ns_per_iter,
    }; 

    if constexpr (DO_STATS)
    {
        let stats = Benchmark_Stats{
            .history_size = history_size,
            .avg_history_size = history_size_running_sum / total_samples,
            .max_history_size = history_size_highwater_mark,
            .min_history_size = history_size_lowwater_mark,

            .delta = required_delta,
            .avg_delta = delta_running_sum / total_samples,
            .max_delta = delta_highwater_mark,
            .min_delta = delta_lowwater_mark,
        };

        return Benchmark_Combined_Results{
            .results = reuslts,
            .stats = stats,
        };
    }
    else
        return reuslts;
}

struct Gather_Samples_Result
{
    std::uint64_t reported_time;
    std::uint64_t total_time;
    std::uint64_t reported_iters;
    std::uint64_t total_iters;

    std::uint64_t block_size;
    std::uint64_t sample_count;
};

template <thunk Fn>
func gather_samples(
    detail::ms optimal_time, detail::ms min_time, detail::ms max_time, 
    std::uint64_t timing_samples[], detail::ms samples_size, Fn fn) -> Gather_Samples_Result
{
    using namespace detail;
    using detail::Benchmark_Stage;
    using u64 = std::uint64_t;

    const ms optimal_ms = optimal_time;
    const ms min_ms = min_time;
    const ms max_ms = max_time;

    const time_point start = detail::clock();
    const time_point last = start;
    const time_point after_min = start + min_ms;

    u64 orientation_iters = 0;

    while(true)
    {
        fn();
        orientation_iters ++;
        if(detail::clock() > after_min)
            break;
    }

    const ns orientation_time = start - detail::clock();
    const u64 block_size = [&]{
        u64 ret = 0;
        if(optimal_ms > orientation_time)
        {
            u64 diff = detail::ns_diff(orientation_time, optimal_ms);
            ret = orientation_iters * diff / (orientation_time.count() * samples_size);
        }

        return max(ret, 1);
    }();

    u64 total_samples = 0;
    const time_point after_max = detail::clock() + max_ms;
    while(true)
    {
        for(u64 i = 0; i < block_size; i++)
            fn();

        const time_point now = detail::clock();
        const u64 current_iter_time = detail::ns_diff(last, now);
        timing_samples[total_samples] = current_iter_time;

        total_samples ++;
        last = now;

        if(total_samples >= optimal_samples || now > after_max)
            break;
    }

    const ns reported_time = after_max - detail::clock();
    const ns total_time = orientation_time + reported_time;
    const u64 reported_iters = total_samples * block_size;
    const u64 total_iters = reported_iters + orientation_iters;

    return Gather_Samples_Result{
        .reported_time = reported_time.count(),
        .total_time = total_time.count(),
        .reported_iters = reported_iters,
        .total_iters = total_iters,

        .block_size = block_size,
        .sample_count = total_samples,
    }
}

template <bool DO_STATS = false, bool DO_HISTORY_CHANGE = true, thunk Fn = int(*)()>
func benchmark_avg(std::uint64_t optimal_time, std::uint64_t min_time, std::uint64_t max_time, std::uint64_t optimal_samples, Fn fn, Benchmark_Constants constants) -> std::conditional_t<DO_STATS, Benchmark_Combined_Results, Benchmark_Results>
{
    using detail::Benchmark_Stage;

    min_time *= 1000'000;
    max_time *= 1000'000;
    optimal_time *= 1000'000;

    let sum_geometrict_series = [](double damping){
        return 1.0 / (1.0 - damping);
    };

    constants.history_hardening_grow = 0.3;
    constants.history_softening_grow = 0.3;
    constants.base_history = 0.9;
    constants.history_max = 1;
    constants.history_min = 0.0001;

    constexpr bool DO_DELTA_CHANGE = true;

    double required_delta = constants.base_delta;
    double history_keep = constants.base_history;
    double history_size = sum_geometrict_series(history_keep);

    double history_size_highwater_mark = 0.0;
    double history_size_lowwater_mark = 999999;
    double history_size_running_sum = 0;
    double delta_highwater_mark = 0.0;
    double delta_lowwater_mark = 999999;
    double delta_running_sum = 0;

    std::uint64_t block_size = 1;
    std::uint64_t total_iters = 0;
    std::uint64_t total_time = 0;
    std::uint64_t accepted_samples = 0;
    std::uint64_t total_samples = 0;
    std::uint64_t last_iters = 0;
    std::uint64_t to_time = min_time;
    Benchmark_Stage stage = Benchmark_Stage::Initial;

    std::uint64_t reported_iters = 0;
    std::uint64_t reported_time = 0;
    double running_average_delta = 0;
    double running_average = 0;
    double prev_time = 0.0;

    let start = detail::clock();
    mut last = start;

    //Algorhitm:
    // (1) Run the given function untill min_time: This will lets us calculate block size 
    //     (ammount of times to run the function in sequence counting as a single sample)
    // 
    // (2) Run block and check its runtime difference against previous run. 
    //     If the difference is under some delta we add the sample and make the delta smaller
    //     If the difference is over delta we ignore the sample and increase delta 
    //     We also add the curren difference into the running_average_delta and damp it by history_keep
    //       so that we have sample to check the next runs against.
    //     We also vary ammount of damping of the running_average_delta (history) which lets us control how
    //       far into the past we are considering results. This is important for functions with big variations in runtime.
    //     
    // (3) We report back the obtained data and statistics

    while(true)
    {
        for(std::uint64_t i = 0; i < block_size; i++)
            fn();

        total_iters += block_size;

        let now = detail::clock();
        const std::uint64_t ellapsed = detail::ns_diff(start, now);
        if(ellapsed > to_time)
        {
            if(ellapsed > max_time)
                break;

            if(stage == Benchmark_Stage::Initial)
            {
                total_time += ellapsed;
                if(optimal_time > ellapsed)
                {
                    block_size = total_iters * (optimal_time - ellapsed) / (ellapsed * optimal_samples);
                    block_size = max(block_size, 1);
                }
                else
                    block_size = 1;

                to_time = 0;
                stage = Benchmark_Stage::Blocked;
            }
            else if(stage == Benchmark_Stage::Blocked)
            {   
                //calculate the current block stats
                const std::uint64_t current_iter_time = detail::ns_diff(last, now);
                const std::uint64_t current_iter_count = total_iters - last_iters;

                //calculate delta squared against the last block
                const double curr_time = cast(double) current_iter_time / current_iter_count;
                const double delta = curr_time - running_average;
                const double delta2 = delta * delta;
                const double curr_time2 = curr_time * curr_time;

                running_average = running_average * history_keep + curr_time * (1 - history_keep);

                //if delta is in required range the block is added else is ignored
                const double relative_average_delta = abs(delta / running_average);
                if(relative_average_delta < required_delta)
                {
                    reported_iters += current_iter_count;
                    reported_time += current_iter_time;
                    accepted_samples++;

                    //harden requirements for delta and increase history size => less blocks will get accepted
                    if constexpr (DO_DELTA_CHANGE)
                        required_delta *= constants.delta_hardening_grow;;
                    if constexpr (DO_HISTORY_CHANGE)
                        history_keep += (constants.history_max - history_keep) * constants.history_hardening_grow;
                }
                else
                {
                    //soften requirements for delta and decrease history size => more blocks will get accepted
                    if constexpr (DO_DELTA_CHANGE)
                        required_delta *= constants.delta_softening_grow;;
                    if constexpr (DO_HISTORY_CHANGE)
                        history_keep += (constants.history_min - history_keep) * constants.history_softening_grow;
                }

                //update constants and track stats
                if constexpr (DO_HISTORY_CHANGE)
                    history_size = sum_geometrict_series(history_keep);

                if constexpr (DO_HISTORY_CHANGE && DO_STATS)
                {
                    history_size_running_sum += history_size;
                    history_size_highwater_mark = max(history_size_highwater_mark, history_size);
                    history_size_lowwater_mark = min(history_size_lowwater_mark, history_size);
                }

                if constexpr (DO_DELTA_CHANGE && DO_STATS)
                {
                    delta_running_sum += required_delta;
                    delta_highwater_mark = max(delta_highwater_mark, required_delta);
                    delta_lowwater_mark = min(delta_lowwater_mark, required_delta);
                }

                prev_time = curr_time;
                total_samples++;
                total_time += current_iter_time;
            }

            last = detail::clock();
            last_iters = total_iters;
        }
    }

    if(reported_iters == 0)
    {
        reported_iters = total_iters;
        reported_time = detail::ns_diff(start, detail::clock());
    }

    const double ns_per_iter = cast(double) reported_iters / cast(double) reported_time;
    const double ms_per_iter = ns_per_iter / 1'000'000;
    const std::uint64_t total_real_time = detail::ns_diff(start, detail::clock());

    let reuslts = Benchmark_Results{

        .total_real_time = total_real_time,  
        .total_time = total_time,  
        .total_iters = total_iters,  

        .block_size = block_size,
        .total_samples = total_samples,
        .accepted_samples = accepted_samples,

        .iters = reported_iters,
        .time = reported_time,

        .iters_per_ms = 1 / ms_per_iter,
        .ms_per_iter = ms_per_iter,

        .iters_per_ns = 1 / ns_per_iter,
        .ns_per_iter = ns_per_iter,
    }; 

    if constexpr (DO_STATS)
    {
        let stats = Benchmark_Stats{
            .history_size = history_size,
            .avg_history_size = history_size_running_sum / total_samples,
            .max_history_size = history_size_highwater_mark,
            .min_history_size = history_size_lowwater_mark,

            .delta = required_delta,
            .avg_delta = delta_running_sum / total_samples,
            .max_delta = delta_highwater_mark,
            .min_delta = delta_lowwater_mark,
        };

        return Benchmark_Combined_Results{
            .results = reuslts,
            .stats = stats,
        };
    }
    else
        return reuslts;
}


#undef func 
#undef let 
#undef mut
#undef cast