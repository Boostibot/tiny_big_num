#pragma once

#include <iostream>
#include <limits>

#include "benchmark.h"

#include "big_int.h"
#include "big_int_ops.h"

#include <cstdlib>
#include <cmath>
#include <chrono>

void print(const char* name, Benchmark_Combined_Results const& combined)
{
    auto reuslt = combined.results;
    auto stats = combined.stats;
    std::cout << name << std::endl;
    std::cout << "total_time: " << reuslt.total_real_time / 1'000'000 << "\n";
    std::cout << "ms_per_iter: " << reuslt.ms_per_iter << "\n";
    std::cout << "time: " << reuslt.time/ 1'000'000  << "\n";
    std::cout << "iters: " << reuslt.iters << "\n";
    std::cout << "block_size: " << reuslt.block_size << "\n";
    std::cout << "total_samples: " << reuslt.total_samples << "\n";
    std::cout << "accepted_samples: " << reuslt.accepted_samples << "\n";
    std::cout << "history_size: " << stats.history_size << "\n";
    std::cout << "\tavg: " << stats.avg_history_size << "\n";
    std::cout << "\tmax: " << stats.max_history_size << "\n";
    std::cout << "\tmin: " << stats.min_history_size << "\n";
    std::cout << "delta: " << stats.delta << "\n";
    std::cout << "\tavg: " << stats.avg_delta << "\n";
    std::cout << "\tmax: " << stats.max_delta << "\n";
    std::cout << "\tmin: " << stats.min_delta << "\n";
    std::cout << std::endl;
}


struct Distribution
{
    double mean = 0;
    double deviation = 0;
    double min = 0;
    double max = 0;
    u64 samples = 0;
};

constexpr double INF = std::numeric_limits<double>::infinity();

template <std::ranges::forward_range R>
func get_dsitribution(R const& range, auto selector)
{
    Distribution dist;
    dist.max = -INF;
    dist.min = INF;
    dist.samples = std::ranges::size(range);
    for(auto const& val : range)
        dist.mean += cast(double) selector(val);

    dist.mean /= dist.samples;
    for(auto const& val : range)
    {
        auto selected = cast(double) selector(val);
        auto diff = selected - dist.mean;
        dist.deviation += diff * diff;

        dist.max = max(selected, dist.max);
        dist.min = min(selected, dist.min);
    }

    dist.deviation /= dist.samples;
    dist.deviation = sqrt(dist.deviation);
    return dist;
}


void println(const char* text)
{
    std::cout << text << "\n";
}

std::ostream& operator <<(std::ostream& stream, Distribution dist)
{
    return stream << dist.mean << " [" << dist.min << ", " << dist.max << "] {" << dist.deviation / dist.mean << "}";
}

template <typename Fn, typename Return, typename... Args>
concept callable_ = std::is_invocable_r_v<Return, Fn, Args...>;

//template <typename Fn, typename Fn_Type>
//concept callable = std::is_invocable_r_v<Return, Fn, Args...>;


template <bool DO_STATS = false, bool DO_HISTORY_CHANGE = true, thunk Fn = int(*)()>
func benchmark_avg(detail::u64 optimal_time_ms, detail::u64 min_time_ms, detail::u64 max_time_ms, detail::u64 optimal_samples, Fn fn, Benchmark_Constants constants) -> std::conditional_t<DO_STATS, Benchmark_Combined_Results, Benchmark_Results>
{
    using namespace detail;

    let sum_geometrict_series = [](f64 damping){
        return 1.0 / (1.0 - damping);
    };

    //constexpr bool DO_HISTORY_CHANGE = true;
    constexpr bool DO_DELTA_CHANGE = true;

    const ms optimal_time = ms(optimal_time_ms); 
    const ms min_time = ms(min_time_ms); 
    const ms max_time = ms(max_time_ms);

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
    Benchmark_Stage stage = Benchmark_Stage::Initial;

    u64 reported_iters = 0;
    u64 reported_time = 0;
    f64 running_average = 0;
    f64 prev_time = 0.0;

    time_point start = detail::clock();
    time_point last = start;
    time_point to_time = start + min_time;
    time_point end_time = start + max_time;

    while(true)
    {
        for(u64 i = 0; i < block_size; i++)
            fn();

        total_iters += block_size;

        let now = detail::clock();
        if(now > to_time)
        {
            if(now > end_time)
                break;

            if(stage == Benchmark_Stage::Initial)
            {
                const u64 ellapsed = ns_diff(start, now);
                total_time += ellapsed;

                if(cast(u64) optimal_time.count() > ellapsed)
                {
                    u64 den = max(ellapsed * optimal_samples, 1);
                    block_size = max(total_iters * (optimal_time.count() - ellapsed) / den, 1);
                }
                else
                    block_size = 1;

                stage = Benchmark_Stage::Blocked;
            }
            else if(stage == Benchmark_Stage::Blocked)
            {   
                const u64 current_iter_time = ns_diff(last, now);
                const u64 current_iter_count = total_iters - last_iters;

                const f64 curr_time = cast(f64) current_iter_time / current_iter_count;
                const f64 delta = curr_time - running_average;

                running_average = running_average * history_damping + curr_time * (1 - history_damping);

                const f64 relative_average_delta = abs(delta / running_average);
                if(relative_average_delta < required_delta)
                {
                    reported_iters += current_iter_count;
                    reported_time += current_iter_time;
                    accepted_samples++;

                    if constexpr (DO_DELTA_CHANGE)
                        required_delta *= constants.delta_hardening_grow;;
                    if constexpr (DO_HISTORY_CHANGE)
                        history_damping += (constants.history_min - history_damping) * constants.history_hardening_grow;
                }
                else
                {
                    if constexpr (DO_DELTA_CHANGE)
                        required_delta *= constants.delta_softening_grow;;
                    if constexpr (DO_HISTORY_CHANGE)
                        history_damping += (constants.history_max - history_damping) * constants.history_softening_grow;
                }

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
        reported_time = ns_diff(start, detail::clock());
    }

    const f64 ns_per_iter = cast(f64) reported_iters / cast(f64) reported_time;
    const f64 ms_per_iter = ns_per_iter / 1'000'000;
    const u64 total_real_time = ns_diff(start, detail::clock());

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


double run_test()
{

    using Big_Int = Big_Int_<u64>;
    constexpr size_t MIN_NUM_SIZE = 10000;
    constexpr size_t MAX_NUM_SIZE = 10000;
    constexpr size_t NUM_COUNT = 20;
    constexpr size_t WARM_ITERS = 100;
    constexpr size_t ITERS = 1000'0;

    let gen_random = [](size_t min_num_size, size_t max_num_size) -> Big_Int {
        size_t size = min(max_num_size, max(min_num_size, rand()));

        Big_Int gened;
        resize(&gened, size);
        for(size_t i = 0; i < size; i++)
        {
            const u64 r = rand();
            gened[i] = r ^ (r >> 16) ^ (r >> 32) ^ (r >> 48);
        }

        return gened;
    };

    Big_Int to = gen_random(MAX_NUM_SIZE, MAX_NUM_SIZE);
    Big_Int nums[NUM_COUNT];
    for(size_t i = 0; i < NUM_COUNT; i++)
        nums[i] = gen_random(MIN_NUM_SIZE, MAX_NUM_SIZE);

    let addition = [&]{
        span<u64> to_s = to;
        for(size_t i = 0; i < NUM_COUNT; i++)
        {
            span<const u64> num_s = nums[i];
            cast(void) add_assign_carry(&to_s, num_s);
        }
    };


    i64 val1 = 1;
    i64 val2 = 3;
    i64 iter = 0;
    let constant = [&]{
        int count = 1000;
        for(int i = 0; i < count; i++)
        {
            val1 = 3*val2 + 2;
            val2 = val1/3 + 1;
            iter++;
        }
    };

    let random = [&]{
        int count = rand() % 1000 * 1000;
        for(int i = 0; i < count; i++)
        {
            val1 = 3*val2 + 2;
            val2 = val1/3 + 1;
            iter++;
        }
    };

    size_t time = 500;
    {
        constexpr size_t sets = 8;
        constexpr size_t repeats = 10;
        Benchmark_Combined_Results stats[sets][repeats];
        const char* names[sets] = {
            "RANDOM_LESS_IFS",
            "RANDOM_LESS_IFS_NO_HISTORY_SIZE",
            "RANDOM_LESS_IFS_NO_HISTORY_CHANGE",
            "RANDOM_LESS_IFS_AVG",
            "CONSTANTS_LESS_IFS",
            "CONSTANTS_LESS_IFS_NO_HISTORY_SIZE",
            "CONSTANTS_LESS_IFS_NO_HISTORY_CHANGE",
            "CONSTANTS_LESS_IFS_AVG",
        };
        double avg = 0;
        double varience = 0;

        for(size_t i = 0; i < repeats; i++)
        {
            stats[0][i] = benchmark<true>(time, time / 8, time, 200, random, Benchmark_Constants{});
            stats[1][i] = benchmark<true, true, false>(time, time / 8, time, 200, random, Benchmark_Constants{});
            stats[2][i] = benchmark<true, false>(time, time / 8, time, 200, random, Benchmark_Constants{});
            stats[3][i] = benchmark_avg<true>(time, time / 8, time, 200, random, Benchmark_Constants{});

            stats[4][i] = benchmark<true>(time, time / 8, time, 200, constant, Benchmark_Constants{});
            stats[5][i] = benchmark<true, true, false>(time, time / 8, time, 200, constant, Benchmark_Constants{});
            stats[6][i] = benchmark<true, false>(time, time / 8, time, 200, constant, Benchmark_Constants{});
            stats[7][i] = benchmark_avg<true>(time, time / 8, time, 200, constant, Benchmark_Constants{});
        }

        for(size_t j = 0; j < sets; j++)
        {
            std::cout << "\n" << names[j] << "\n";
            auto stats_span = span<Benchmark_Combined_Results>{stats[j], repeats};

            std::cout << "ns_per_iter:\t" << 
                get_dsitribution(stats_span, [](auto res){ return res.results.ns_per_iter; }) << "\n";
            std::cout << "iters:\t" 
                << get_dsitribution(stats_span, [](auto res){ return res.results.iters; }) << "\n";
            std::cout << "avg_delta:\t" 
                << get_dsitribution(stats_span, [](auto res){ return res.stats.avg_delta; }) << "\n";
            std::cout << "min_delta:\t" 
                << get_dsitribution(stats_span, [](auto res){ return res.stats.min_delta; }) << "\n";
            std::cout << "max_delta:\t" 
                << get_dsitribution(stats_span, [](auto res){ return res.stats.max_delta; }) << "\n";
            std::cout << "accepted/rejected:\t" 
                << get_dsitribution(stats_span, [](auto res){ return res.results.accepted_samples / cast(double) res.results.total_samples; }) << "\n";
        }
    }


    //std::cout << "\nADD:" << std::endl;
    //benchmark(2000, addition);
    //benchmark(2000, addition);
    //benchmark(2000, addition);
    //benchmark(2000, addition);

    return 1.0;
}