#pragma once


#include <iostream>
#include <cctype>
#include <random>
#include <limits>
#include <initializer_list>
#include <memory_resource>
#include <vector>
#include <array>
#include <algorithm>

//#define USE_CUSTOM_LIB

#ifdef USE_CUSTOM_LIB
#include "jot/stack.h"
#include "jot/allocator_arena.h"

#undef move
#undef forward

#endif // USE_CUSTOM_LIB

#include "benchmark.h"
#include "big_int_ops.h"


namespace test
{

    using Unsigned_Max = u64;
    using std::size;
    using std::data;

    #ifdef USE_CUSTOM_LIB
        template <typename T>
        using Vector = jot::Stack<T, 16>;

        using Memory_Resource = jot::Allocator_Resource;

        using Monotonic_Buffer_Resource = jot::Unbound_Arena_Resource;

        template<typename T>
        func make_sized_vector(size_t size, Memory_Resource* resource) -> Vector<T>
        {
            using namespace jot;
            Vector<T> vec = Poly_Allocator{std::move(resource)};
            if(size != 0)
                force(resize(&vec, cast(isize) size));
            return vec;
        }

        runtime_proc release_memory(Monotonic_Buffer_Resource* resource) noexcept -> void 
        {
            resource->deallocate_all();
        }

        template<typename T>
        func push_back(Vector<T>* vec, T pushed)
        {
            using namespace jot;
            force(push(vec, move(pushed)));
        }

        template<typename T>
        func resize(Vector<T>* vec, size_t size) -> void
        {
            jot::force(jot::resize(vec, size));
        }

    #else
        template <typename T>
        using Vector = std::pmr::vector<T>;

        using Memory_Resource = std::pmr::memory_resource;

        using Monotonic_Buffer_Resource = std::pmr::monotonic_buffer_resource;

        template<typename T>
        func make_sized_vector(size_t size, Memory_Resource* resource) -> Vector<T>
        {
            Vector<T> vec(std::move(resource));
            if(size != 0)
                vec.resize(size);
            return vec;
        }


        runtime_proc release_memory(Monotonic_Buffer_Resource* resource) noexcept -> void
        {
            resource->release();
        }

        template<typename T>
        func push_back(Vector<T>* vec, T pushed)
        {
            vec->push_back(move(pushed));
        }

        template<typename T>
        func resize(Vector<T>* vec, size_t size)
        {
            vec->resize(size);
        }
    #endif

    template<typename T>
    func to_slice(Vector<T> const& vec) -> Slice<const T>
    {
        return Slice<const T>{data(vec), cast(usize) size(vec)};
    }

    template<typename T>
    func to_slice(Vector<T>* vec) -> Slice<T>
    {
        return Slice<T>{std::data(*vec), cast(usize) size(*vec)};
    }



    template<typename T>
    func make_vector_of_slice(Slice<const T> const& slice, Memory_Resource* resource) -> Vector<T>
    {
        mut vec = make_sized_vector<T>(slice.size, resource);
        mut vec_slice = to_slice(&vec);
        copy_slice<T>(&vec_slice, slice, Iter_Direction::ANY);

        return vec;
    }

    template<typename T>
    func make_vector_of_digits(Unsigned_Max val, Memory_Resource* resource) -> Vector<T>
    {
        constexpr size_t count = digits_to_represent<T, Unsigned_Max>();
        Vector<T> vec = make_sized_vector<T>(count, resource);
        Slice<T> slice = to_slice(&vec);
        Slice<T> converted_slice = from_number<T>(&slice, val);
        resize(&vec, converted_slice.size);
        return vec;
    }

    template <typename T>
    struct Padded_Vector
    {
        Vector<T> vector;
        size_t prefix_size;
        size_t content_size;
        size_t postfix_size;
    };

    template <typename T>
    func make_padded_vector_of_digits(Unsigned_Max val, Memory_Resource* resource, size_t pref_size = 1, size_t post_size = 1) -> Padded_Vector<T>
    {
        constexpr size_t count = digits_to_represent<T, Unsigned_Max>();
        size_t total_size = pref_size + count + post_size;

        Padded_Vector<T> padded;

        padded.vector = make_sized_vector<T>(total_size, resource);
        Slice<T> whole = to_slice(&padded.vector);
        Slice<T> digits = slice_size(whole, pref_size, count);
        Slice<T> content = from_number<T>(&digits, val);

        padded.prefix_size = pref_size;
        padded.content_size = content.size;
        padded.postfix_size = post_size;

        return padded;
    }

    template <typename T>
    func content(Padded_Vector<T>* padded) -> Slice<T>
    {
        Slice<T> whole = to_slice(&padded->vector);
        return slice_size(whole, padded->prefix_size, padded->content_size);
    }

    template <typename T>
    func content(Padded_Vector<T> const& padded) -> Slice<const T>
    {
        Slice<const T> whole = to_slice(padded.vector);
        return slice_size(whole, padded.prefix_size, padded.content_size);
    }

    using Random_Generator = std::mt19937;
    struct Batch_Op_Result
    {
        Unsigned_Max value;
        Unsigned_Max overflow;

        bool constexpr operator ==(Batch_Op_Result const&) const noexcept = default;
    };

    template <typename T>
    constexpr size_t MAX_TYPE_SIZE_FRACTION = sizeof(Unsigned_Max) / sizeof(T);

    constexpr size_t RELEASE_MEMORY_EVERY = 10;

    template <typename T>
    struct Uniform_Exponential_Distribution
    {
        double shift;
        double quotient; //the range of values for which the distribution should yield the same probability
        std::uniform_real_distribution<> uniform_distribution;

        Uniform_Exponential_Distribution(T from, T to, T quotient)
        {
            *this = Uniform_Exponential_Distribution(cast(double) from, cast(double) to, cast(double) quotient);
        }

        Uniform_Exponential_Distribution(double from, double to, double quotient)
        {
            assert(quotient != 1 && quotient > 0 && "quotient must be bigger than 0 and not 1");
            assert(to >= from && "from .. to must be a valid range");
            double shift = 1 - from; //in case from is 0 or negative

            assert(to + shift >= 1);
            double log_quot = log(quotient);
            double from_log = log(from + shift) / log_quot;
            double to_log = log(to + shift) / log_quot;

            this->shift = shift;
            this->quotient = quotient;
            this->uniform_distribution = std::uniform_real_distribution<>{from_log, to_log};
        }

        template<class Generator>
        T operator()(Generator& gen) const 
        {
            double val = pow(this->quotient, this->uniform_distribution(gen)) - this->shift;
            if constexpr (std::is_floating_point_v<T>)
                return cast(T) val;
            else
                return cast(T) round(val);
        }

    };

    template <typename T, typename ReturnT = Unsigned_Max>
    runtime_func make_exponential_distribution(ReturnT min = 0, ReturnT max = std::numeric_limits<ReturnT>::max()) -> Uniform_Exponential_Distribution<ReturnT>
    {
        double quot = pow(2, BIT_SIZE<T>);
        return {cast(double) min, cast(double) max, quot};
    }

    runtime_func generate_random_optims(Random_Generator* gen) -> Optim_Info
    {
        let exp_dist = make_exponential_distribution<u8, size_t>(0, 100);
        //std::uniform_int_distribution<> dist{0, 100};
        std::uniform_int_distribution<> uniform_bool{0, 1};
        let random_bool = [&]() -> bool { return cast(bool) uniform_bool(*gen); };
        let random_num = [&]() -> size_t { return cast(size_t) exp_dist(*gen); };

        Optim_Info optims;
        optims.mul_shift = random_bool();
        optims.mul_consts = random_bool();
        optims.mul_half_bits = random_bool();

        optims.div_shift = random_bool();
        optims.div_consts = random_bool();
        optims.rem_optims = random_bool();

        optims.max_reusion_depth = random_num();
        optims.mul_quadratic_both_below_size = random_num();
        optims.mul_quadratic_single_below_size = random_num();
        optims.trivial_pow_below_power = random_num();

        return optims;
    }

    u64 ipow_squares_native(u64 x, u64 n)
    {
        if (x <= 1) 
            return x;

        u64 y = 1;
        for (; n != 0; n >>= 1, x *= x)
            if (n & 1)
                y *= x;

        return y;
    }

    u64 ipow_squares(u64 x, u64 n)
    {
        if (x <= 1) 
            return x;
        if(n == 0)  
            return 1;

        u64 i = 0;
        u64 y = 1;
        size_t max_pos = find_last_set_bit(n);
        for(; i <= max_pos; i++)
        {
            i64 bit = get_bit<u64>(n, i);
            if(bit)
                y *= x;
            x *= x;
        }
        return y;
    }

    u64 ipow_trivial(u64 x, u64 n)
    {
        u64 res = 1;
        while(n > 0)
        {
            res *= x;
            n--;
        }

        return res;
    }

    u64 iroot_shifting(u64 x, u64 n)
    {
        if(n == 0)
            return 1;
        if (x <= 1) 
            return x;

        u64 r = 1;
        u64 s = ((find_last_set_bit(x) / n) * n);

        while(true)
        {
            if(s < n) //fast
                break;

            s -= n; //fast
            r <<= 1; //O(n)
            u64 power = ipow_squares(r | 1, n); //slow
            u64 bit = power <= (x >> s);
            r |= bit;
        }

        return r;
    }

    u64 iroot_newton(u64 of_value, u64 power)
    {
        const u64 x = of_value; 
        const u64 n = power;

        if(x == 0)
            return 1;
        if (x <= 1) 
            return x;

        //u64 u = x;
        //u64 s = x+1;

        // We want to find such r that:
        //   r^n = x
        // this can be rewritten as: 
        //   n*log(r) = log(x)
        // factoring r gives:
        //   r = e^(log(x) / n)
        //   r = 2^(log2(x) / 2)

        // we want to find an initial estimate above such r so we need to increase the expression:
        //  we will perform this by parts:
        //  log2(x) <= [log2(x)] + 1
        //  log2(x)/n <= upper[ ([log2(x)] + 1) / n ] = [ ([log2(x)] + 1 + n - 1) / n ]
        //            = [ ([log2(x)] + n) / n ]
        // 
        //  so the upper estimate for r is:
        //      2^[ ([log2(x)] + n) / n ]

        const u64 log2x = find_last_set_bit(x);
    
        const u64 uppe_initial_estimate = cast(u64) 1 << ((log2x + n) / n);
        u64 r = uppe_initial_estimate;
        u64 prev_r = -1;

        while(true)
        {
            // the newton update
            u64 new_upper_r = (n-1) * r + x / ipow_squares(r, n-1);

            prev_r = r;
            r = new_upper_r / n;

            //we monotonically decrease 
            // => prev_r should always be bigger than r
            // => if not we break and return the last value that 
            //     satisfied the condition (prev_r)
            if(r >= prev_r)
                break;
        }

        return prev_r;
    }

    u64 ilog_heyley(u64 of_value, u64 base)
    {
        const u64 x = of_value;
        const u64 b = base;

        //b^y = x
        //y = log_b(x)

        if(x <= 1)
            return 0;
        if(base <= 1)
            return 0;

        const u64 log2x = find_last_set_bit(x);
        const u64 log2b = find_last_set_bit(b);

        assert(log2x != 0);
        assert(log2b != 0);

        const u64 lower_initial_estimate = log2x / (log2b + 1);
        const u64 higher_initial_estimate = (log2x + log2b) / (log2b);
        u64 y = higher_initial_estimate;
        u64 prev_y = y;

        while(true)
        {
            u64 value_from_curr_approximate = ipow_squares(b, y);
            u64 c_x = value_from_curr_approximate;

            prev_y = y;
            if(c_x <= x)
                break;

            u64 new_y = y - 2*(c_x - x) / (x + c_x);


            //if we didnt move at all we are one above the answer (is equivalent to the commented out code
            if(new_y == y)
                return new_y - 1;

            //if(new_y == y) 
                //new_y --;
            
            y = new_y;
        }

        return prev_y;

    }

    u64 ilog_bsearch(u64 of_value, u64 base)
    {
        const u64 x = of_value;
        const u64 b = base;
    
        if(x <= 1)
            return 0;
        if(base <= 1)
            return 0;

        const u64 log2x = find_last_set_bit(x);
        const u64 log2b = find_last_set_bit(b);

        assert(log2x != 0);
        assert(log2b != 0);

        const u64 lower_initial_estimate = log2x / (log2b + 1);
        const u64 higher_initial_estimate = (log2x + log2b) / (log2b);

        u64 curr_lo = lower_initial_estimate;
        u64 curr_hi = higher_initial_estimate;

        u64 y = (curr_lo + curr_hi) / 2;

        while(true)
        {
            if(curr_hi - curr_lo <= 1)
                break;

            u64 value_from_curr_approximate = ipow_squares(b, y);
            u64 c_x = value_from_curr_approximate;

            if(c_x > x)
                curr_hi = y;
            else if(c_x < x)
                curr_lo = y;
            else
            {
                curr_lo = y;
                break;
            }

            y = (curr_lo + curr_hi) / 2;
        }

        return curr_lo;
    }

    runtime_proc test_misc(Random_Generator* generator, size_t random_runs)
    {
        using Max = Unsigned_Max;

        assert(find_last_set_bit<Max>(0) == 0);
        assert(find_last_set_bit<Max>(1) == 0);
        assert(find_last_set_bit<Max>(2) == 1);
        assert(find_last_set_bit<Max>(3) == 1);
        assert(find_last_set_bit<Max>(4) == 2);
        assert(find_last_set_bit<Max>(7) == 2);
        assert(find_last_set_bit<Max>(8) == 3);
        assert(find_last_set_bit<Max>(76) == 6);
        assert(find_last_set_bit<Max>(513) == 9);
        assert(find_last_set_bit<Max>(641) == 9);
        assert(find_last_set_bit<Max>(1024) == 10);
        assert(find_last_set_bit<u64>(0xf4f1f63db9b7e800) == 63);
        assert(find_last_set_bit<Max>(0b010100) == 4);
        assert(find_last_set_bit<Max>(0b010000) == 4);

        assert(find_first_set_bit<Max>(0) == 0);
        assert(find_first_set_bit<Max>(0b000001) == 0);
        assert(find_first_set_bit<Max>(0b010100) == 2);
        assert(find_first_set_bit<Max>(0b010101) == 0);
        assert(find_first_set_bit<Max>(0b010000) == 4);
        assert(find_first_set_bit<Max>(0b010100) == 2);
        assert(find_first_set_bit<Max>(0b010010) == 1);

        assert(pop_count<Max>(0b0000'0000) == 0);
        assert(pop_count<Max>(0b0000'0001) == 1);
        assert(pop_count<Max>(0b1000'0000) == 1);
        assert(pop_count<Max>(0b1010'0000) == 2);
        assert(pop_count<Max>(0b1010'0011) == 4);
        assert(pop_count<Max>(0b0111'1111) == 7);
        assert(pop_count<Max>(cast(u32) -1) == 32);
        assert(pop_count<Max>(cast(Max) -1) == 64);

        assert(high_mask<u8>() == 0xF0);
        assert(high_mask<u8>(2) == 0b1111'1100);
        assert(high_mask<u8>(5) == 0b1110'0000);
        assert(high_mask<u16>() == 0xFF00);
        assert(high_mask<u32>() == 0xFFFF0000);
        assert(high_mask<u64>() == 0xFFFFFFFF00000000);

        assert(low_mask<u8>() == 0x0F);
        assert(low_mask<u8>(2) == 0b0000'0011);
        assert(low_mask<u8>(5) == 0b0001'1111);
        assert(low_mask<u16>() == 0x00FF);
        assert(low_mask<u32>() == 0x0000FFFF);
        assert(low_mask<u64>() == 0x00000000FFFFFFFF);

        assert(low_bits<u32>(0b1010'1100'0001, 1) == 1);
        assert(low_bits<u32>(0b1010'1100'0001, 2) == 1);
        assert(low_bits<u32>(0b1010'1100'0001, 4) == 1);
        assert(low_bits<u32>(0b1010'1100'0001, 6) == 1);
        assert(low_bits<u32>(0b1010'1100'0001, 8) == 0b1100'0001);
        assert(low_bits<u32>(1, 16) == 1);
        assert(low_bits<u32>(1, 7) == 1);
        assert(low_bits<u32>(1, 8) == 1);


        assert(set_bit(0b0001001, 1, 1) == 0b0001011);
        assert(set_bit(0b0001001, 0, 1) == 0b0001001);
        assert(set_bit(0b0001001, 0, 0) == 0b0001000);
        assert(set_bit(0b0001001, 3, 0) == 0b0000001);
        assert(get_bit(0b0001001, 0) == 1);
        assert(get_bit(0b0001001, 1) == 0);
        assert(get_bit(0b0001001, 2) == 0);
        assert(get_bit(0b0001001, 3) == 1);


        assert(ipow_squares(3,3) == 27);

        assert(iroot_shifting(8, 3) == 2);
        assert(iroot_shifting(9, 3) == 2);
        assert(iroot_shifting(27, 3) == 3);
        assert(iroot_shifting(28, 3) == 3);
        assert(iroot_shifting(29, 3) == 3);
        assert(iroot_shifting(34, 5) == 2);
        assert(iroot_shifting(15625, 3) == 25);
        assert(iroot_shifting(0, 0) == 1);
        assert(iroot_shifting(1, 0) == 1);
        assert(iroot_shifting(4096, 6) == 4);
        assert(iroot_shifting(18446744073709551614, 17) == 13);

        assert(iroot_newton(8, 3) == 2);
        assert(iroot_newton(9, 3) == 2);
        assert(iroot_newton(27, 3) == 3);
        assert(iroot_newton(28, 3) == 3);
        assert(iroot_newton(29, 3) == 3);
        assert(iroot_newton(34, 5) == 2);
        assert(iroot_newton(15625, 3) == 25);
        assert(iroot_newton(0, 0) == 1);
        assert(iroot_newton(1, 0) == 1);
        assert(iroot_newton(4096, 6) == 4);


        assert(ilog_heyley(25, 2) == 4);
        assert(ilog_heyley(255, 2) == 7);
        assert(ilog_heyley(256, 2) == 8);
        assert(ilog_heyley(4, 4) == 1);
        assert(ilog_heyley(9, 3) == 2);
        assert(ilog_heyley(26, 3) == 2);
        assert(ilog_heyley(27, 3) == 3);
        assert(ilog_heyley(64, 3) == 3);
        assert(ilog_heyley(81, 3) == 4);
        assert(ilog_heyley(16, 4) == 2);
        assert(ilog_heyley(17, 4) == 2);
        assert(ilog_heyley(27, 4) == 2);
        assert(ilog_heyley(63, 4) == 2);
        assert(ilog_heyley(64, 4) == 3);


        assert(ilog_bsearch(25, 2) == 4);
        assert(ilog_bsearch(255, 2) == 7);
        assert(ilog_bsearch(256, 2) == 8);
        assert(ilog_bsearch(4, 4) == 1);
        assert(ilog_bsearch(9, 3) == 2);
        assert(ilog_bsearch(26, 3) == 2);
        assert(ilog_bsearch(27, 3) == 3);
        assert(ilog_bsearch(64, 3) == 3);
        assert(ilog_bsearch(81, 3) == 4);
        assert(ilog_bsearch(16, 4) == 2);
        assert(ilog_bsearch(17, 4) == 2);
        assert(ilog_bsearch(27, 4) == 2);
        assert(ilog_bsearch(63, 4) == 2);
        assert(ilog_bsearch(64, 4) == 3);


        constexpr u64 max = std::numeric_limits<u64>::max()/2;
        constexpr u64 min = 2;
        let val_dist = make_exponential_distribution<u64>(min, max);

        assert(cast(u64) 1 << 42 == 4398046511104);
        assert(ilog_heyley(4398046511104, 2) == 42);

        for(size_t i = 0; i < random_runs; i++)
        {
            u64 val = cast(u64) val_dist(*generator);
            let power_dist = make_exponential_distribution<u64>(min, val);
            u64 power = cast(u64) power_dist(*generator);

            u64 rooted = iroot_shifting(val, power);
            u64 powed_low = ipow_squares(rooted, power);
            u64 back_down = iroot_shifting(powed_low, power);
            if(rooted > 1 && powed_low < 43980465111) //so that no overflow
            {
                u64 loged_heyley = ilog_heyley(powed_low, rooted);
                u64 loged_bsearch = ilog_heyley(powed_low, rooted);
                assert(loged_heyley == power);
                assert(loged_bsearch == power);
            }

            if(power < 100)
            {
                u64 ref_powed_low = ipow_trivial(rooted, power);
                assert(ref_powed_low == powed_low);
            }
            else
            {
                //most of the time (maybe always) the diff is exactly zero 
                // yet we still want to be sure this test wont break randomly
                double epsilon = 0.00000000001;
                double ref_ipow = pow(rooted, power);
                double powed_low_d = cast(double) powed_low;

                double diff = abs(ref_ipow - powed_low_d);
                double frac = diff / cast(double) val; 
                //we divide by val because the only way errors can happen is double conversion on
                // large numbers => we need to only consider large numbers and their 'relative' error
                assert(frac < epsilon);
            }
            assert(powed_low <= val);
            assert(back_down == rooted);


            u64 rooted_again = iroot_shifting(rooted, power);
            u64 rooted_again_newton = iroot_newton(rooted, power);

            assert(rooted_again == rooted_again_newton);
        }
    }

    template <typename T>
    runtime_proc test_add_overflow(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;
        
        runtime_proc test_add_overflow_batch = [&](Max left_, Max right_, Max carry_in, Res expected, bool check_overflow = true) -> bool {
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Padded pad_right = make_padded_vector_of_digits<T>(right_, resource, 1);

            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);
            Slice pad_right_s = to_slice<T>(&pad_right.vector);

            let loc = Op_Location::OUT_OF_PLACE;

            let res1 = ::add_overflow_batch<T>(&out_s,  to_slice(left), to_slice(right), loc, cast(T) carry_in);
            let res2 = ::add_overflow_batch<T>(&pad_left_s, content(pad_left), to_slice(right), loc, cast(T) carry_in); //test self assign
            let res3 = ::add_overflow_batch<T>(&pad_right_s, to_slice(left), content(pad_right), loc, cast(T) carry_in); 

            Slice res1_s = res1.slice;
            Slice res2_s = res2.slice;
            Slice res3_s = res3.slice;

            if(check_overflow == false)
            {
                //@TODO: make not shitty
                //add the overflow as extra digit and use the result to obtain the converted number
                res1_s = trim(out_s, res1.slice.size + 1);
                res2_s = trim(pad_left_s, res2.slice.size + 1);
                res3_s = trim(pad_right_s, res3.slice.size + 1);

                *last(&res1_s) = res1.overflow;
                *last(&res2_s) = res2.overflow;
                *last(&res3_s) = res3.overflow;
            }

            let num1 = unwrap(to_number(res1_s));
            let num2 = unwrap(to_number(res2_s));
            let num3 = unwrap(to_number(res3_s));

            bool outputs_match = 
                num1 == expected.value && 
                num2 == expected.value && 
                num3 == expected.value;

            bool overflows_match =
                res1.overflow == res2.overflow && 
                res1.overflow == res3.overflow;

            bool sizes_match = 
                res1.slice.size == res2.slice.size && 
                res1.slice.size == res3.slice.size;

            if(check_overflow)
                overflows_match = overflows_match && res1.overflow == expected.overflow;

            return outputs_match && overflows_match && sizes_match;
        };
    
        runtime_proc auto_test_add = [&](Max left, Max right, Max carry_in = 0) -> bool {
            Max highest_one = cast(Max) 1 << (BIT_SIZE<Max> - 1);
            Max adjusted_left = left & ~highest_one;
            Max adjusted_right = right & ~highest_one;
            Max adjusted_carry = cast(Max) low_bits(cast(T) carry_in);
            Max expected = add_no_overflow(adjusted_left, adjusted_right, adjusted_carry);

            return test_add_overflow_batch(adjusted_left, adjusted_right, adjusted_carry, Res{expected, 0}, false);
        };

        runtime_proc test_single_add = [&](Max left, Max single){
            if(left == 0 || single == 0)
                return true;

            Max highest_one = cast(Max) 1 << (BIT_SIZE<Max> - 1);
            Max adjusted_left = left & ~highest_one;
            T adjusted_single = cast(T) (single & ~highest_one);

            Max expected = adjusted_left + adjusted_single;

            mut left1 = make_vector_of_digits<T>(adjusted_left, resource); 
            mut left2 = make_vector_of_digits<T>(adjusted_left, resource);

            Slice left1_s = to_slice(&left1);
            Slice left2_s = to_slice(&left2);

            let res2 = add_overflow_batch<T>(&left1_s, left1_s, adjusted_single, Op_Location::IN_PLACE);
            let res1 = add_overflow_batch<T>(&left2_s, left2_s, Slice{&adjusted_single, 1}, Op_Location::IN_PLACE);

            bool slices = is_equal<T>(res1.slice, res2.slice);
            bool overflows = res1.overflow == res2.overflow;
            return slices && overflows;
        };

        if constexpr(std::is_same_v<T, u8>)
        {
            assert(test_add_overflow_batch(0, 0, 0, Res{0, 0}));
            assert(test_add_overflow_batch(0, 1, 0, Res{1, 0}));
            assert(test_add_overflow_batch(1, 1, 0, Res{2, 0}));

            assert(test_add_overflow_batch(0x0025, 0x00FF, 0, Res{0x24, 1}));
            assert(test_add_overflow_batch(0x00FF, 0x0025, 0, Res{0x24, 1}));
            assert(test_add_overflow_batch(0xFFFFFFFF, 0xFFFFFFFF, 0, Res{0xFFFFFFFE, 1}));

            assert(test_add_overflow_batch(0x0025, 0x01FF, 0, Res{0x224, 0}));
            assert(test_add_overflow_batch(0x01FF, 0x0025, 0, Res{0x224, 0}));
            assert(test_add_overflow_batch(0x01FF, 0, 0, Res{0x01FF, 0}));
            assert(test_add_overflow_batch(0x01FF, 1, 0, Res{0x0200, 0}));
            assert(test_add_overflow_batch(0x01FF, 1, 1, Res{0x0201, 0}));
            assert(test_add_overflow_batch(1, 0x01FF, 0, Res{0x0200, 0}));

            assert(test_add_overflow_batch(0x1225, 0xEFEF, 0, Res{0x0214, 1}));


            release_memory(&local_buffer);
        }

        assert(auto_test_add(0xFF, 0xFF));
        assert(auto_test_add(0xFFFF'FFFF'FFFF'FFFF, 0xFFFF'FFFF'FFFF'FFFF));
        assert(auto_test_add(56775421, 32101454));
        assert(auto_test_add(64344, 45));
        assert(auto_test_add(984643131, 81766));
        assert(auto_test_add(18748341324, 13875));
        assert(auto_test_add(23446634457799, 454406999));
        assert(auto_test_add(74984193, 982387980));

        assert(test_single_add(130, cast(T) 196));
        assert(test_single_add(74, cast(T) 173));

        let distribution = make_exponential_distribution<T>();
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Max left = distribution(*generator);
            Max right = distribution(*generator);
            Max carry = distribution(*generator);
            assert(auto_test_add(left, right, carry));
            assert(test_single_add(left, cast(T) right));
        };
    }

    
    template <typename T>
    runtime_proc test_sub_overflow(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        runtime_proc test_sub_overflow_batch = [&](Max left_, Max right_, Max carry_in, Res expected, bool check_overflow = true) -> bool {
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Padded pad_right = make_padded_vector_of_digits<T>(right_, resource, 1);

            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);

            let loc = Op_Location::OUT_OF_PLACE;

            let res1 = ::sub_overflow_batch<T>(&out_s,  to_slice(left), to_slice(right), loc, cast(T) carry_in);
            let res2 = ::sub_overflow_batch<T>(&pad_left_s, content(pad_left), to_slice(right), loc, cast(T) carry_in); //test self assign

            Slice res1_s = res1.slice;
            Slice res2_s = res2.slice;

            if(check_overflow == false)
                assert(left_ >= right_ && "must not overflow when auto testing");

            let num1 = unwrap(to_number(res1_s));
            let num2 = unwrap(to_number(res2_s));

            bool outputs_match =    num1 == expected.value && num2 == expected.value;
            bool overflows_match =  res1.overflow == res2.overflow;
            bool sizes_match =      res1.slice.size == res2.slice.size;

            if(check_overflow)
                overflows_match = overflows_match && res1.overflow == expected.overflow;

            return outputs_match && overflows_match && sizes_match;
        };

        runtime_proc auto_test_sub = [&](Max left, Max right, Max carry_in = 0) -> bool {
            if(left < right)
                std::swap(left, right);

            Max highest_one = cast(Max) 1 << (BIT_SIZE<Max> - 1);
            Max adjusted_right = right & ~highest_one;
            Max adjusted_carry = carry_in > low_mask<Max>();

            Max adjusted_left = left;
            if(left < adjusted_carry + adjusted_right)
                adjusted_left = add_no_overflow(left, adjusted_carry);

            Max expected = adjusted_left - adjusted_right - adjusted_carry;
            return test_sub_overflow_batch(adjusted_left, adjusted_right, adjusted_carry, Res{expected, 0});
        };


        if constexpr(std::is_same_v<T, u8>)
        {
            assert(test_sub_overflow_batch(0, 0, 0, Res{0, 0}));
            assert(test_sub_overflow_batch(1, 0, 0, Res{1, 0}));
            assert(test_sub_overflow_batch(0, 1, 0, Res{0xFF, 1}));
            assert(test_sub_overflow_batch(1, 1, 0, Res{0, 0}));
            assert(test_sub_overflow_batch(2, 1, 0, Res{1, 0}));

            assert(test_sub_overflow_batch(0xFF0000, 0x000001, 0, Res{0xFEFFFF, 0}));
            assert(test_sub_overflow_batch(0x000000, 0x000001, 0, Res{0xFF, 1}));
            assert(test_sub_overflow_batch(0xFF0000, 0x0000FF, 0, Res{0xFEFF01, 0}));

            assert(test_sub_overflow_batch(0xFFFF, 0xFF, 0, Res{0xFF00, 0}));
            assert(test_sub_overflow_batch(0xFFFE, 0xFF, 0, Res{0xFEFF, 0}));
            assert(test_sub_overflow_batch(0xFE, 0xFF, 0, Res{0xFF, 1}));
            assert(test_sub_overflow_batch(0x0, 0xFFFFFFFF, 0, Res{0x1, 1}));
            assert(test_sub_overflow_batch(0x1, 0xFFFFFFFF, 0, Res{0x2, 1}));
            assert(test_sub_overflow_batch(0xFFFFFFFF, 0x1, 0, Res{0xFFFFFFFE, 0}));

            release_memory(&local_buffer);
        }

        assert(auto_test_sub(0xFF, 0xFF));
        assert(auto_test_sub(0xFFFF'FFFF'FFFF'FFFF, 0xFFFF'FFFF'FFFF'FFFF));
        assert(auto_test_sub(56775421, 32101454));
        assert(auto_test_sub(64344, 45));
        assert(auto_test_sub(984643131, 81766));
        assert(auto_test_sub(18748341324, 13875));
        assert(auto_test_sub(23446634457799, 454406999));
        assert(auto_test_sub(74984193, 982387980));
        assert(auto_test_sub(0x04b941000665a200, 0x7d57e, 0xc87afb6c));

        let distribution = make_exponential_distribution<T>();
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Max left = distribution(*generator);
            Max right = distribution(*generator);
            Max carry = distribution(*generator);
            assert(auto_test_sub(left, right, carry));
        };
    }

    template <typename T>
    runtime_proc test_complement_overflow(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        //we dont test oveflow of this op 
        runtime_proc complement_overflow_batch = [&](Max left_, Max carry_in = 1) -> Max{
            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Slice left_s = content(&pad_left);
            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);

            let loc = Op_Location::OUT_OF_PLACE;

            let res1 = ::complement_overflow_batch<T>(&out_s, left_s, cast(T) carry_in);
            let res2 = ::complement_overflow_batch<T>(&pad_left_s, left_s, cast(T) carry_in); //test self assign

            let num1 = unwrap(to_number(res1.slice));
            let num2 = unwrap(to_number(res2.slice));

            assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
            assert(num1 == num2 && "self assign must work correctly");

            return num1;
        };


        runtime_proc auto_test_complement = [&](Max left, Max carry_in = 1) -> bool{
            Max complemented = ~left + carry_in;
            Max cropped_complemented = complemented;

            size_t bits = 0;
            if(left == 0)
                bits = 0;
            else if(BIT_SIZE<T> >= BIT_SIZE<Max>)
                bits = BIT_SIZE<Max>;
            else
            {
                for(Max temp_left = left; temp_left > 0;)
                {
                    temp_left >>= BIT_SIZE<T>;
                    bits += BIT_SIZE<T>;
                }
            }
        
            if(bits < BIT_SIZE<Max>)
                cropped_complemented = low_bits(complemented, bits);

            return complement_overflow_batch(left, carry_in) == cropped_complemented;
        };
    
        if constexpr(std::is_same_v<T, u8>)
        {
            assert(complement_overflow_batch(0) == 0);
            assert(complement_overflow_batch(1) == cast(u8)(-1));
            assert(complement_overflow_batch(5) == cast(u8)(-5));
            assert(complement_overflow_batch(0xFF) == cast(u8)(-0xFF));
            assert(complement_overflow_batch(0xABCD) == cast(u16)(-0xABCD));
            assert(complement_overflow_batch(0xFFFF0000) == cast(u32)(-cast(i64)0xFFFF0000));
        }

        let exponential = make_exponential_distribution<T, Max>();
        let uniform = std::uniform_int_distribution<>(0, 1);
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Max left = exponential(*generator);
            Max carry = uniform(*generator);
            assert(auto_test_complement(left, carry));
        };
    }

    template <typename T>
    runtime_proc test_shift_overflow(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        runtime_proc shift_both = [&](Max left_, Max right, Max carry_in, bool up_down, Iter_Direction direction) -> Batch_Op_Result{
            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Slice left_s = content(&pad_left);
            Slice out_s = to_slice<T>(&out);
            Slice whole = to_slice<T>(&pad_left.vector);

            Slice pad_left_s;
            if(direction == Iter_Direction::FORWARD)
                pad_left_s = whole;
            else
                pad_left_s = slice(whole, 1);

            Batch_Overflow<T> res1;
            Batch_Overflow<T> res2;

            if(up_down)
            {
                res1 = ::shift_up_overflow_batch<T>(&out_s, left_s, right, direction, cast(T) carry_in);
                res2 = ::shift_up_overflow_batch<T>(&pad_left_s, left_s, right, direction, cast(T) carry_in);
            }
            else
            {
                res1 = ::shift_down_overflow_batch<T>(&out_s, left_s, right, direction, cast(T) carry_in);
                res2 = ::shift_down_overflow_batch<T>(&pad_left_s, left_s, right, direction, cast(T) carry_in);
            }

            let num1 = unwrap(to_number(res1.slice));
            let num2 = unwrap(to_number(res2.slice));

            assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
            assert(res1.overflow == res2.overflow);
            assert(num1 == num2 && "self assign must work correctly");

            return Batch_Op_Result{num1, res1.overflow};
        };

        runtime_proc shift_up_overflow_batch = [&](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
            let forward = shift_both(left, right, carry_in, true, Iter_Direction::FORWARD);
            let backward = shift_both(left, right, carry_in, true, Iter_Direction::BACKWARD);

            assert(forward.value == backward.value);
            assert(forward.overflow == backward.overflow);

            return forward;
        };

        runtime_proc shift_down_overflow_batch = [&](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
            let forward = shift_both(left, right, carry_in, false, Iter_Direction::FORWARD);
            let backward = shift_both(left, right, carry_in, false, Iter_Direction::BACKWARD);

            assert(forward.value == backward.value);
            assert(forward.overflow == backward.overflow);

            return forward;
        };


        if constexpr(std::is_same_v<T, u8>)
        {
            assert((shift_up_overflow_batch(0, 0, 0) == Res{0, 0}));
            assert((shift_up_overflow_batch(0, 1, 0) == Res{0, 0}));
            assert((shift_up_overflow_batch(0, 1, 0xFF) == Res{0, 0xFF})); //should be consistent with multiply

            assert((shift_up_overflow_batch(1, 1) == Res{2, 0}));
            assert((shift_up_overflow_batch(1, 4) == Res{16, 0}));
            assert((shift_up_overflow_batch(1, 4, 0b1010'0000) == Res{0b11010, 0}));
            assert((shift_up_overflow_batch(11, 4) == Res{176, 0}));
            assert((shift_up_overflow_batch(0xAB, 0x4) == Res{0xB0, 0xA}));
            assert((shift_up_overflow_batch(0xAB, 0x4, 0xD0) == Res{0xBD, 0xA}));
            assert((shift_up_overflow_batch(0xAB, 0x4, 0x0D) == Res{0xB0, 0xA}));
            assert((shift_up_overflow_batch(0b0101'0101'0111'0001, 7) == Res{0b1011'1000'1000'0000, 0b010'1010}));
            assert((shift_up_overflow_batch(0b0101'0101'0111'0001, 0) == Res{0b0101'0101'0111'0001, 0}));
            assert((shift_up_overflow_batch(0b0101'0101'0111'0001, 7, 0b11001100) == Res{0b1011'1000'1110'0110, 0b010'1010}));

            release_memory(&local_buffer);

            assert((shift_down_overflow_batch(0, 0, 0) == Res{0, 0}));
            assert((shift_down_overflow_batch(0, 1, 0) == Res{0, 0}));
            assert((shift_down_overflow_batch(0, 1, 0xFF) == Res{0, 0xFF}));
            assert((shift_down_overflow_batch(1, 1) == Res{0, 0x80}));
            assert((shift_down_overflow_batch(0b10, 1) == Res{0b1, 0}));
            assert((shift_down_overflow_batch(0b110000, 3) == Res{0b110, 0}));
            assert((shift_down_overflow_batch(0b0101'0011, 1) == Res{0b101001, 0x80}));
            assert((shift_down_overflow_batch(0b0101'0011, 4) == Res{0b101, 0b0011'0000}));
            assert((shift_down_overflow_batch(0b0101'0011, 4, 0b0000'1111) == Res{0b1111'0101, 0b0011'0000}));
            assert((shift_down_overflow_batch(0b0101'0011, 0) == Res{0b0101'0011, 0}));
            assert((shift_down_overflow_batch(0b0101'0011, 7) == Res{0b0, 0b10100110}));
            assert((shift_down_overflow_batch(0b0101'0011, 7, 0b1100'0111) == Res{0b1000'1110, 0b10100110}));
            assert((shift_down_overflow_batch(0xFFEEDD, 4) == Res{0x0FFEED, 0xD0}));
            assert((shift_down_overflow_batch(0xFFEEDD, 4, 0x0A) == Res{0xAFFEED, 0xD0}));
            assert((shift_down_overflow_batch(0b0101010101110001, 5) == Res{0b0000001010101011, 0b10001000}));
            assert((shift_down_overflow_batch(0b0101010101110001, 5, 0b1010'1111) == Res{0b0111101010101011, 0b10001000}));

            release_memory(&local_buffer);
        }
    }

    struct Adjusted
    {
        Unsigned_Max left;
        Unsigned_Max right;
    };

    runtime_func adjust_to_not_mul_overflow(Unsigned_Max left, Unsigned_Max right) -> Adjusted
    {
        using Max = Unsigned_Max;
        // (by doing find_last_set_bit) and then shifting it so that they will not overflow

        Max adjusted_left = left;
        Max adjusted_right = right;

        size_t log_left = find_last_set_bit(left) + 1; 
        size_t log_right = find_last_set_bit(right) + 1; 
        //find_last_set_bit returns position of highest set bit but we want "under which power of two is the whole number" 
        // => we add 1

        size_t combined_log = log_left + log_right;

        //if can overflow we shift down so that its not possible
        if(combined_log > BIT_SIZE<Max>)
        {
            size_t diff = combined_log - BIT_SIZE<Max>;
            size_t shift = (diff + 1) / 2;

            adjusted_left >>= shift;
            adjusted_right >>= shift;
        }

        assert(find_last_set_bit(adjusted_left) + find_last_set_bit(adjusted_right) + 2 <= BIT_SIZE<Max>);
        return {adjusted_left, adjusted_right};
    }

    template <typename T>
    runtime_proc test_mul_overflow(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        runtime_proc test_mul_overflow_batch = [&](Max left_, Max right, Max carry_in, Res expected, Optim_Info optims, bool check_overflow = true) -> bool {
            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Slice left_s = content(&pad_left);
            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);

            let res1 = ::mul_overflow_batch<T>(&out_s, left_s, cast(T) right, optims, cast(T) carry_in);
            let res2 = ::mul_overflow_batch<T>(&pad_left_s, left_s, cast(T) right, optims, cast(T) carry_in);

            Slice res1_s = res1.slice;
            Slice res2_s = res2.slice;

            if(check_overflow == false)
            {
                //add the overflow as extra digit and use the result to obtain the converted number
                res1_s = trim(out_s, res1.slice.size + 1);
                res2_s = trim(pad_left_s, res2.slice.size + 1);

                *last(&res1_s) = res1.overflow;
                *last(&res2_s) = res2.overflow;
            }

            let num1 = unwrap(to_number(res1_s));
            let num2 = unwrap(to_number(res2_s));

            bool outputs_match = 
                num1 == expected.value && 
                num2 == expected.value; 

            bool overflows_match = res1.overflow == res2.overflow; 
            bool sizes_match = res1.slice.size == res2.slice.size; 

            if(check_overflow)
                overflows_match = overflows_match && res1.overflow == expected.overflow;

            return outputs_match && overflows_match && sizes_match;
        };

        runtime_proc auto_test_mul = [&](Max left, Max right, Optim_Info optims) -> bool
        {
            let adjusted = adjust_to_not_mul_overflow(left, right);
            return test_mul_overflow_batch(adjusted.left, adjusted.right, 0, Res{adjusted.left * adjusted.right, 0}, optims, false);
        };

        if constexpr(std::is_same_v<T, u8>)
        {
            for(size_t i = 0; i < controlled_runs; i++)
            {
                let optims = generate_random_optims(generator);
                assert(test_mul_overflow_batch(0, 0, 0, Res{0, 0}, optims));
                assert(test_mul_overflow_batch(1, 0, 0, Res{0, 0}, optims));
                assert(test_mul_overflow_batch(0, 1, 0, Res{0, 0}, optims));
                assert(test_mul_overflow_batch(1, 1, 0, Res{1, 0}, optims));
                assert(test_mul_overflow_batch(0, 1, 0xFF, Res{0, 0xFF}, optims));
                assert(test_mul_overflow_batch(25, 1, 0, Res{25, 0}, optims));
                assert(test_mul_overflow_batch(25, 5, 0, Res{125, 0}, optims));
                assert(test_mul_overflow_batch(0xFF, 0x10, 0, Res{0xF0, 0xF}, optims));

                assert((mul_overflow<T>(0xFF, 0x10, 0xC) == Overflow<T>{0xFC, 0xF}));
                assert(test_mul_overflow_batch(0xFFEE, 0x80, 0, Res{0xF700, 0x7F}, optims));
                assert(test_mul_overflow_batch(0xFFEE, 0x80, 0, Res{0xF700, 0x7F}, optims));
                assert(test_mul_overflow_batch(0xA1, 0x11, 0, Res{0xB1, 0xA}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 1, 0, Res{0xABCDEF00, 0}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 0, 0, Res{0, 0}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 0xAB, 0, Res{0xC28EA500, 0x72}, optims));
                assert(test_mul_overflow_batch(13745, 137, 0, Res{0xBBB9, 0x1C}, optims));


                release_memory(&local_buffer);
            }
        }

        assert(auto_test_mul(0x14156f, 0x3c, Optim_Info{}));

        let exponential = make_exponential_distribution<T, Max>();
        let uniform = std::uniform_int_distribution<long long>(0, FULL_MASK<T> >> 1);
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Max left = exponential(*generator);
            Max right = uniform(*generator);
            let optims = generate_random_optims(generator);
            assert(auto_test_mul(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_div_overflow_low(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        assert((div_overflow_low<T>(0x12, 0x0A, 0) == Overflow<T>{1, 0x08}));

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        runtime_proc test_div_overflow_low_batch = [&](Max left_, Max right, Max carry_in, Res expected, Optim_Info optims) -> bool {
            constexpr size_t padding_before = 2;

            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Slice left_s = content(&pad_left);
            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = slice(to_slice<T>(&pad_left.vector), 1);

            let res1 = ::div_overflow_low_batch<T>(&out_s, left_s, cast(T) right, optims, cast(T) carry_in);
            let res2 = ::div_overflow_low_batch<T>(&pad_left_s, left_s, cast(T) right, optims, cast(T) carry_in);

            let num1 = unwrap(to_number(res1.slice));
            let num2 = unwrap(to_number(res2.slice));

            bool outputs_match = 
                num1 == expected.value && 
                num2 == expected.value; 

            bool overflows_match = 
                res1.overflow == expected.overflow && 
                res2.overflow == expected.overflow; 

            bool sizes_match = res1.slice.size == res2.slice.size; 

            return outputs_match && overflows_match && sizes_match;
        };
    
        runtime_proc auto_test_div = [&](Max left, Max right, Optim_Info optims) -> bool {
            if(right == 0)
                return true;
        
            Max adjusted_right = low_bits(cast(T) right);
            Max res = left / adjusted_right;
            Max rem = left % adjusted_right;

            //carry in is basically added as an extra digit to the divided number
            // its entire purpose is linking multiple div blocks if desired (ie we store number chunked for some reason)
            // we dont test this yet because the setup would basically exactly mirror of the operation
            // (because its dependent on the template type)

            return test_div_overflow_low_batch(left, adjusted_right, 0, Res{res, rem}, optims);
        };

        if constexpr(std::is_same_v<T, u8>)
        {
            for(size_t i = 0; i < controlled_runs; i++)
            {
                let optims = generate_random_optims(generator);
                assert(test_div_overflow_low_batch(0, 1, 0, Res{0, 0}, optims));
                assert(test_div_overflow_low_batch(1, 1, 0, Res{1, 0}, optims));
                assert(test_div_overflow_low_batch(25, 1, 0, Res{25, 0}, optims));
                assert(test_div_overflow_low_batch(25, 5, 0, Res{5, 0}, optims));
                assert(test_div_overflow_low_batch(0xFF, 0x08, 0, Res{31, 7}, optims));
                assert(test_div_overflow_low_batch(0xA1, 0x0F, 0, Res{10, 11}, optims));
                assert(test_div_overflow_low_batch(0xA1C0, 0x08, 0x1, Res{0x3438, 0}, optims));
                assert(test_div_overflow_low_batch(0xABCDEF00, 1, 0, Res{0xABCDEF00, 0}, optims));
                assert(test_div_overflow_low_batch(0x99AA, 0x08, 0x9, Res{0x3335, 2}, optims));
                assert(test_div_overflow_low_batch(0xABCDEF00, 0x07, 0xC, Res{0xcf668fdb, 0x3}, optims));
                //for carry in the result is obtained as follows:
                // left:  0x0ABCDEF00 
                // carry: 0xC00000000
                // eff:   0xCABCDEF00 
                // eff / right = 0x1CF668FDB => 0xCF668FDB (removal of the added - carried digit)
                // eff % right = 0x3
                assert(test_div_overflow_low_batch(0xABCDEF00, 0x07, 0, Res{411771428, 4}, optims));
                assert(test_div_overflow_low_batch(0xABCDEF00, 0x0F, 0, Res{192160000, 0}, optims));
                assert(test_div_overflow_low_batch(13745, 0xB, 0, Res{1249, 6}, optims));

                release_memory(&local_buffer);
            }
        }

        let exponential = make_exponential_distribution<T, Max>();
        let uniform = std::uniform_int_distribution<long long>(0, low_mask<T>());
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Max left = exponential(*generator);
            Max right = uniform(*generator);
            let optims = generate_random_optims(generator);
            assert(auto_test_div(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_mul_quadratic(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        runtime_proc test_fused_mul_add = [&](Max left_, Max right_, Max coef, Max expected, Optim_Info optims = Optim_Info{}) -> bool{
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);
            Vector output = make_sized_vector<T>(size(left), resource);

            Slice left_s = to_slice(&left);
            Slice right_s = to_slice(&right);
            Slice output_s = to_slice(&output);

            let res_out = ::fused_mul_add_overflow_batch<T>(&output_s, left_s, right_s, cast(T) coef, optims, Op_Location::OUT_OF_PLACE);
            let res_in = ::fused_mul_add_overflow_batch<T>(&left_s, left_s, right_s, cast(T) coef, optims, Op_Location::IN_PLACE);
        
            push_back(&output, res_out.overflow);
            push_back(&left, res_in.overflow);

            left_s = to_slice(&left);
            output_s = to_slice(&output);

            let num_out = unwrap(to_number(left_s));
            let num_in = unwrap(to_number(output_s));

            return expected == num_out && expected == num_in; 
        };

        runtime_proc auto_test_fused_mul_add = [&](Max left, Max right, Max coef, Optim_Info optims = Optim_Info{}) -> bool{
            return test_fused_mul_add(left, right, coef, left + coef*right, optims);
        };

        runtime_proc test_mul_quadratic = [&](Max left_, Max right_, Max expected, Optim_Info optims = Optim_Info{}) -> bool{
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);

            size_t to_size = required_mul_to_size(size(left), size(right));
            size_t aux_size1 = required_mul_quadratic_auxiliary_size(size(left), size(right));
            size_t aux_size2 = required_mul_karatsuba_auxiliary_size(size(left), size(right));

            size_t aux_size = max(aux_size1, aux_size2)*10; //to allow all karatsuba recursions

            Vector res = make_sized_vector<T>(to_size, resource);
            Vector aux = make_sized_vector<T>(aux_size, resource);

            Slice res_s = to_slice(&res);
            Slice aux_s = to_slice(&aux);
            CSlice left_s = to_slice(left);
            CSlice right_s = to_slice(right);

            ::mul_quadratic<T>(&res_s, &aux_s, left_s, right_s, optims);
            let normal = unwrap(to_number(res_s));

            ::mul_quadratic_fused<T>(&res_s, left_s, right_s, optims);
            let fused = unwrap(to_number(res_s));

            ::mul_karatsuba<T>(&res_s, &aux_s, left_s, right_s, optims);
            let karatsuba = unwrap(to_number(res_s));

            return expected == normal && expected == fused && expected == karatsuba;
        };

        runtime_proc auto_test_mul_quadratic = [&](Max left, Max right, Optim_Info optims) -> bool
        {
            let adjusted = adjust_to_not_mul_overflow(left, right);
            return test_mul_quadratic(adjusted.left, adjusted.right, adjusted.left * adjusted.right, optims);
        };

        Optim_Info shift = {true, true, false};

        assert(test_fused_mul_add(1, 1, 1, 2));
        assert(test_fused_mul_add(1, 1, 1, 2, shift));
        assert(test_fused_mul_add(1, 1, 5, 6));
        assert(test_fused_mul_add(1, 1, 5, 6, shift));
        assert(test_fused_mul_add(10, 1, 5, 15, shift));
        assert(auto_test_fused_mul_add(0xFFFF, 1, 0x80, shift));
        assert(auto_test_fused_mul_add(0x1000000000, 0x574ce80, 0x10, shift));
        assert(auto_test_fused_mul_add(0x0da8824af2592c00, 0x0da8824af2592c00, 0x2, shift));

        assert(test_mul_quadratic(0, 0, 0));
        assert(test_mul_quadratic(1, 0, 0));
        assert(test_mul_quadratic(0, 1, 0));
        assert(test_mul_quadratic(1, 1, 1));
        assert(test_mul_quadratic(25, 1, 25));
        assert(test_mul_quadratic(25, 5, 125));
        assert(test_mul_quadratic(0xFF, 0x10, 0xFF0));
        assert(test_mul_quadratic(0xA1, 0x11, 0xAB1));
        assert(test_mul_quadratic(0xABCDEF00, 1, 0xABCDEF00));
        assert(test_mul_quadratic(0xABCDEF00, 0, 0));
        assert(test_mul_quadratic(13745, 137, 0x1CBBB9));
        assert(test_mul_quadratic(0xFFFF, 0x10000, 0xFFFF0000));
        assert(test_mul_quadratic(0xABCD, 0x10001, 0xABCDABCD));
        assert(test_mul_quadratic(0xAB, 0x1254, 0xAB*0x1254));
        assert(test_mul_quadratic(0x1000, 0x1200, 0x1000*0x1200));
        assert(test_mul_quadratic(0xABCD, 0x54, 0xABCD*0x54));
        assert(test_mul_quadratic(0xABCD, 0x1254, 0xABCD*0x1254));
        assert(test_mul_quadratic(0xABCDEF124, 0x1254, 0xC4CDA61BA7D0));
        assert(test_mul_quadratic(0x1254, 0xABCDEF124, 0xC4CDA61BA7D0));
        release_memory(&local_buffer);

        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x38, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0xFFFF, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x1010, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x7838, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x27838, shift));
        assert(auto_test_mul_quadratic(60136640465, 262154, shift));

        let distribution = make_exponential_distribution<T>();
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Optim_Info optims = generate_random_optims(generator);
            Max left = distribution(*generator);
            Max right = distribution(*generator);

            assert(auto_test_mul_quadratic(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_div_bit_by_bit(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        enum Will_Fail
        {
            FAIL,
            OKAY
        };

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        runtime_proc test_div_bit_by_bit = [&](Max num_, Max den_, Max ex_quotient, Max ex_remainder, Will_Fail expected_fail = OKAY, Optim_Info optims = Optim_Info{}) -> bool{
            Vector num = make_vector_of_digits<T>(num_, resource); 
            Vector den = make_vector_of_digits<T>(den_, resource);
            size_t quo_size = required_div_quotient_size(size(num), size(den));
            size_t rem_size = required_div_remainder_size(size(num), size(den));

            Vector quo = make_sized_vector<T>(quo_size, resource); 
            Vector rem = make_sized_vector<T>(rem_size, resource);

            Slice quo_s = to_slice(&quo);
            Slice rem_s = to_slice(&rem);
            Slice num_s = to_slice(&num);
            Slice den_s = to_slice(&den);

            let div_res = ::div_bit_by_bit<T>(&quo_s, &rem_s, num_s, den_s, optims);
            if(expected_fail == FAIL)
                return div_res.has == false;

            let unwrapped = unwrap(div_res);
            let quotient = unwrap(to_number(unwrapped.quotient));
            let remainder = unwrap(to_number(unwrapped.remainder));
            return ex_quotient == quotient && ex_remainder == remainder && div_res.has;
        };


        runtime_proc auto_test_div_bit_by_bit = [&](Max left, Max right, Optim_Info optims = Optim_Info{}) -> bool{
            if(right == 0)
                return test_div_bit_by_bit(left, right, 0, 0, FAIL);

            Max quotient = left / right;
            Max remainder = left % right;
            return test_div_bit_by_bit(left, right, quotient, remainder, OKAY, optims);
        };

        assert(test_div_bit_by_bit(0, 0, 0, 0, FAIL));
        assert(test_div_bit_by_bit(1, 0, 0, 0, FAIL));
        assert(test_div_bit_by_bit(0, 1, 0, 0, OKAY));
        assert(test_div_bit_by_bit(1, 1, 1, 0));
        assert(test_div_bit_by_bit(25, 1, 25, 0));
        assert(test_div_bit_by_bit(1, 25, 0, 1));
        assert(test_div_bit_by_bit(25, 26, 0, 25));
        assert(test_div_bit_by_bit(25, 5, 5, 0));
        assert(test_div_bit_by_bit(28, 5, 5, 3));
        assert(test_div_bit_by_bit(0xFF, 0x10, 0xF, 0xF));
        assert(test_div_bit_by_bit(183, 27, 6, 21));
        assert(test_div_bit_by_bit(0xABCDEF00, 1, 0xABCDEF00, 0));
        assert(test_div_bit_by_bit(0xABCDEF00, 0xABCDEF00, 1, 0));
        assert(test_div_bit_by_bit(0xABCDEF00, 0xABCDEF10, 0, 0xABCDEF00));
        assert(test_div_bit_by_bit(0xABCDEF00, 0, 0, 0, FAIL));
        assert(test_div_bit_by_bit(0xAB, 0xABCD, 0, 0xAB));
        assert(test_div_bit_by_bit(0xABCE, 0xABCD, 1, 1));
        release_memory(&local_buffer);

        assert(auto_test_div_bit_by_bit(13745, 6435));
        assert(auto_test_div_bit_by_bit(4643477, 4796));
        assert(auto_test_div_bit_by_bit(7893, 01047));
        assert(auto_test_div_bit_by_bit(65536, 256));
        assert(auto_test_div_bit_by_bit(256, 65536));
        assert(auto_test_div_bit_by_bit(256*413541, 256));
        assert(auto_test_div_bit_by_bit(959464934, 1111));
        assert(auto_test_div_bit_by_bit(27417318639, 57432413));
        assert(auto_test_div_bit_by_bit(6994935563924682752, 67877739413342216));

        let distribution = make_exponential_distribution<T>();
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Optim_Info optims = generate_random_optims(generator);
            Unsigned_Max left = distribution(*generator);
            Unsigned_Max right = distribution(*generator);
            assert(auto_test_div_bit_by_bit(left, right, optims));
        };
    }

    template <typename T>
    void println(Slice<T> slice)
    {
        std::cout << "[";
        if(slice.size > 0)
        {
            std::cout << cast(Unsigned_Max) slice[0];

            for(size_t i = 1; i < slice.size; i++)
                std::cout << ", " << cast(Unsigned_Max) slice[i];
        }
        std::cout << "]\n";
    }

    template <typename Num, typename Rep, typename Conv_Fn>
    func to_base(Slice<const Num> num, Num base, Conv_Fn conv, Optim_Info optims, Memory_Resource* resource) -> Vector<Rep> {
        Vector<Num> temp = make_vector_of_slice<Num>(num, resource);
        Vector<Rep> out = make_sized_vector<Rep>(required_size_to_base<Num>(num.size, base), resource);

        Slice<Num> temp_s = to_slice(&temp);
        Slice<Rep> out_s = to_slice(&out);

        Slice<Rep> converted = (to_base<Num, Rep>(&out_s, &temp_s, num, base, conv, optims));
        assert(is_striped_representation(converted));

        resize(&out, converted.size);
        return out;
    };

    template <typename Num, typename Rep, typename Conv_Fn>
    func from_base(Slice<const Rep> rep, Num base, Conv_Fn conv, Optim_Info optims, Memory_Resource* resource) -> Vector<Num> {
        size_t required_size = required_size_from_base<Num>(rep.size, base);
        Vector<Num> out = make_sized_vector<Num>(required_size + 1, resource); //@TODO: why is there the +1

        Slice<Num> out_s = to_slice(&out);
        Slice<Num> converted = from_base<Num, Rep>(&out_s, rep, base, conv, optims);

        assert(is_striped_number(converted));
        resize(&out, converted.size);
        return out;
    };

    func stringlen(const char* str) -> size_t
    {
        if(str == nullptr)
            return 0;

        size_t size = 0;
        while(str[size] != '\0')
            size ++;

        return size;
    }

    template <typename Num = Unsigned_Max>
    func to_base(Slice<const Num> num, Num base, Memory_Resource* resource) -> Vector<char> {
        assert(base > 2);
        assert(base <= 36);

        constexpr char val_to_char_table[36] = {
            '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z'
        };
    
        let conversion = [&](Num value) -> char {
            assert(value < 36);
            return val_to_char_table[value];
        };

        return to_base<Num, char>(num, base, conversion, Optim_Info{}, resource);
    }

    template <typename Num = Unsigned_Max>
    func from_base(const char* str, Num base, Memory_Resource* resource) -> Trivial_Maybe<Vector<Num>> {
        assert(base > 2);
        assert(base <= 36);

        proc make_table = [](){
            std::array<Num, 255> table;

            for(size_t i = 0; i < table.size(); i++)
                table[i] = cast(Num) -1;

            table['0'] = 0;
            table['1'] = 1;
            table['2'] = 2;
            table['3'] = 3;
            table['4'] = 4;
            table['5'] = 5;
            table['6'] = 6;
            table['7'] = 7;
            table['8'] = 8;
            table['9'] = 9;
            table['A'] = 10;
            table['a'] = 10;
            table['B'] = 11;
            table['b'] = 11;
            table['C'] = 12;
            table['c'] = 12;
            table['D'] = 13;
            table['d'] = 13;
            table['E'] = 14;
            table['e'] = 14;
            table['F'] = 15;
            table['f'] = 15;
            table['G'] = 16;
            table['g'] = 16;
            table['H'] = 17;
            table['h'] = 17;
            table['I'] = 18;
            table['i'] = 18;
            table['J'] = 19;
            table['j'] = 19;
            table['K'] = 20;
            table['k'] = 20;
            table['L'] = 21;
            table['l'] = 21;
            table['M'] = 22;
            table['m'] = 22;
            table['N'] = 23;
            table['n'] = 23;
            table['O'] = 24;
            table['o'] = 24;
            table['P'] = 25;
            table['p'] = 25;
            table['Q'] = 26;
            table['q'] = 26;
            table['R'] = 27;
            table['r'] = 27;
            table['S'] = 28;
            table['s'] = 28;
            table['T'] = 29;
            table['t'] = 29;
            table['U'] = 30;
            table['u'] = 30;
            table['V'] = 31;
            table['v'] = 32;
            table['W'] = 32;
            table['w'] = 32;
            table['X'] = 33;
            table['x'] = 33;
            table['Y'] = 34;
            table['y'] = 34;
            table['Z'] = 35;
            table['z'] = 35;

            return table;
        };

        constexpr std::array<Num, 255> char_to_val_table = make_table();

        bool ok = true;
        let conversion = [&](char value) -> Num {
            Num val = char_to_val_table[value];
            if(val == -1)
            {
                ok = false;
                return 0;
            }

            return val;
        };

        Slice<const char> rep{str, stringlen(str)};
        Vector<Num> res = from_base<Num, char>(rep, base, conversion, Optim_Info{}, resource);

        if(ok)
            return wrap(res);
        else
            return {};
    }

    template <typename Num, typename Rep>
    func native_to_base(Unsigned_Max num, Num base, Memory_Resource* resource) -> Vector<Rep> 
    {
        using Max = Unsigned_Max;
        const size_t num_digits = digits_to_represent<Num, Max>();
        const size_t max_size = required_size_to_base<Num>(num_digits, base);

        Vector<Rep> converted = make_sized_vector<Rep>(max_size, resource);
        size_t size = 0;
        for(Max curr_val = num; curr_val != 0; curr_val /= base, size ++)
        {
            Max rem = curr_val % base;
            converted[size] = cast(Rep) rem;
        }

        resize(&converted, size);
        Slice<Rep> converted_s = to_slice(&converted);
        reverse(&converted_s);

        return converted;
    };

    template <typename Num, typename Rep>
    func native_from_base(Slice<const Rep> rep, Num base) -> Unsigned_Max 
    {
        using Max = Unsigned_Max;
        assert(is_striped_representation(rep));

        Max value = 0;
        for(size_t i = 0; i < rep.size; i++)
        {
            Max muled = value * base;
            Max added = muled + rep[i];
            assert(muled >= value);
            assert(added >= muled);

            value = added;
        }

        return value;
    };

    template <typename Num, typename Rep>
    runtime_proc test_to_base(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs)
    {
        static_assert(sizeof(Num) >= sizeof(Rep));
        using Max = Unsigned_Max;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        proc id_to_conversion = [](Num num) -> Rep {return cast(Rep) num;};
        proc id_from_conversion = [](Rep num) -> Num {return cast(Num) num;};

        runtime_proc test_to_base = [&](Slice<const Num> num, Num base, Slice<const Rep> expected_rep, Optim_Info optims) -> bool 
        {
            Vector<Rep> converted = to_base<Num, Rep>(num, base, id_to_conversion, optims, resource);
        
            assert(is_striped_representation(expected_rep));

            return is_equal<Rep>(to_slice(converted), expected_rep);
        };

        runtime_proc manual_test_to_base = [&](Max num, Num base, std::initializer_list<Rep> expected_rep, Optim_Info optims = Optim_Info{}) -> bool 
        {
            Vector<Num> num_ = make_vector_of_digits<Num>(num, resource);
            Slice<const Rep> expected_rep_s = {std::data(expected_rep), std::size(expected_rep)};
            return test_to_base(to_slice(num_), base, expected_rep_s, optims);
        };

        runtime_proc auto_test_to_base = [&](Max num, Num base, Optim_Info optims) -> bool 
        {
            Vector<Num> num_ = make_vector_of_digits<Num>(num, resource);
            Vector<Rep> expected_rep = native_to_base<Num, Rep>(num, base, resource);
            return test_to_base(to_slice(num_), base, to_slice(expected_rep), optims);
        };

        runtime_proc test_from_base = [&](Slice<const Rep> rep, Num base, Slice<const Num> expected_num, Optim_Info optims) -> bool 
        {
            Vector<Num> num = from_base<Num, Rep>(rep, base, id_from_conversion, optims, resource);
        
            assert(is_striped_representation(expected_num));
            assert(is_striped_representation(to_slice(num)));

            return is_equal<Num>(to_slice(num), expected_num);
        };

        runtime_proc manual_test_from_base = [&](std::initializer_list<Rep> rep, Num base, Max expected_num, Optim_Info optims = Optim_Info{}) -> bool 
        {
            Vector<Num> expected_num_ = make_vector_of_digits<Num>(expected_num, resource);
            Slice<const Rep> rep_s = {std::data(rep), std::size(rep)};
            return test_from_base(rep_s, base, to_slice(expected_num_), optims);
        };

        runtime_proc auto_test_from_base = [&](Slice<const Rep> rep, Num base, Optim_Info optims) -> bool 
        {
            Max num = native_from_base<Num, Rep>(rep, base);
            Vector<Num> expected_num = make_vector_of_digits<Num>(num, resource);
            return test_from_base(rep, base, to_slice(expected_num), optims);
        };

        runtime_proc test_to_and_fro = [&](Max value, Num base, Optim_Info optims) -> bool
        {
            Vector<Num> initial = make_vector_of_digits<Num>(value, resource);
            let rep = to_base<Num, Rep>(to_slice(initial), base, id_to_conversion, optims, resource);
            let num = from_base<Num, Rep>(to_slice(rep), base, id_from_conversion, optims, resource);
            Slice<const Num> num_s = to_slice(num);
            let obtained_value = unwrap(to_number(num_s));

            bool match = is_equal<Num>(to_slice(num), to_slice(initial));
            return match;
        };

        assert((manual_test_to_base(1, 2, {1})));
        assert((manual_test_to_base(2, 2, {1,0})));
        assert((manual_test_to_base(4, 2, {1,0,0})));
        assert((manual_test_to_base(5, 2, {1,0,1})));
        assert((manual_test_to_base(10, 2, {1,0,1,0})));

        assert((manual_test_to_base(0, 9, {})));
        assert((manual_test_to_base(1, 9, {1})));
        assert((manual_test_to_base(2, 9, {2})));
        assert((manual_test_to_base(65453, 9, {1,0,8,7,0,5})));
        assert((manual_test_to_base(844354, 9, {1,5,2,6,2,1,1})));
        assert((manual_test_to_base(844354, 10, {8,4,4,3,5,4})));
        release_memory(&local_buffer);

        assert((manual_test_from_base({1}, 2, 1)));
        assert((manual_test_from_base({1,0}, 2, 2)));
        assert((manual_test_from_base({1,0,0}, 2, 4)));
        assert((manual_test_from_base({1,0,1}, 2, 5)));
        assert((manual_test_from_base({1,0,1,0}, 2, 10)));

        assert((manual_test_from_base({}, 9, 0)));
        assert((manual_test_from_base({1}, 9, 1)));
        assert((manual_test_from_base({2}, 9, 2)));
        assert((manual_test_from_base({1,0,8,7,0,5}, 9, 65453)));
        assert((manual_test_from_base({1,5,2,6,2,1,1}, 9, 844354)));
        assert((manual_test_from_base({8,4,4,3,5,4}, 10, 844354)));

        assert((test_to_and_fro(10, 2, Optim_Info{})));
        assert((test_to_and_fro(59180470887883776, 3, Optim_Info{})));

        let num_dist = make_exponential_distribution<Num>();
        constexpr size_t max_1 = std::numeric_limits<Rep>::max() >> HALF_BIT_SIZE<Rep>;
        constexpr size_t max_2 = std::numeric_limits<Num>::max() >> HALF_BIT_SIZE<Num>;

        let base_dist = make_exponential_distribution<Num, Num>(2, cast(Num) min(max_1, max_2));
    
        for(size_t i = 0; i < random_runs; i++)
        {
            if(i % RELEASE_MEMORY_EVERY == 0)
                release_memory(&local_buffer);

            Optim_Info optims = generate_random_optims(generator);
            Unsigned_Max num = num_dist(*generator);
            Num base = base_dist(*generator);

            assert(auto_test_to_base(num, base, optims));
            assert(test_to_and_fro(num, base, optims));
        };
    }

    template <typename T>
    runtime_proc test_pow_by_squaring(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using Max = Unsigned_Max;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Monotonic_Buffer_Resource local_buffer{upstream};
        Memory_Resource* resource = &local_buffer;

        runtime_proc pow_ = [&](CSlice num, size_t pow, Optim_Info optims, bool do_trivial) -> Vector {
            size_t required_size = required_pow_to_size(num, pow);
            size_t aux_size = required_pow_by_squaring_auxiliary_size(num, pow)*100; //for karatsuba
            Vector out = make_sized_vector<T>(required_size, resource);
            Vector aux = make_sized_vector<T>(aux_size, resource);

            Slice out_s = to_slice(&out);
            Slice aux_s = to_slice(&aux);

            Slice powed;
            if(do_trivial)
                powed = ::trivial_pow<T>(&out_s, &aux_s, num, pow, optims);
            else
                powed = ::pow_by_squaring<T>(&out_s, &aux_s, num, pow, optims);

            resize(&out, powed.size);
            return out;
        };

        runtime_proc pow = [&](Max num, size_t pow, Optim_Info optims, bool do_trivial) -> Max {
            Vector num_ = make_vector_of_digits<T>(num, resource);
            Vector powed = pow_(to_slice(num_), pow, optims, do_trivial);
            CSlice powed_s = to_slice(powed);

            Max res = unwrap(to_number(powed_s));
            return res;
        };

        runtime_proc pow_by_squaring = [&](Max num, Max power, Optim_Info optims = Optim_Info{}) -> Max {
            return pow(num, power, optims, false);
        };

        runtime_proc trivial_pow = [&](Max num, Max power, Optim_Info optims = Optim_Info{}) -> Max {
            return pow(num, power, optims, true);
        };

        runtime_proc test_pow = [&](CSlice num, size_t pow, CSlice expected, Optim_Info optims) -> bool {
            Vector powed = pow_(num, pow, optims, false);
            CSlice powed_s = to_slice(powed);

            assert(is_striped_number(powed_s));
            assert(is_striped_number(expected));
        
            return is_equal<T>(powed_s, expected);
        };


        assert(trivial_pow(0, 0) == 1);
        assert(trivial_pow(0, 1) == 0);
        assert(trivial_pow(0, 54) == 0);
        assert(trivial_pow(0, 5554135) == 0);
        assert(trivial_pow(2, 0) == 1);
        assert(trivial_pow(1048576, 0) == 1);
        assert(trivial_pow(5465313, 1) == 5465313);
        assert(trivial_pow(5465313, 2) == 5465313ull*5465313ull);
        assert(trivial_pow(5465, 3) == 5465ull*5465ull*5465ull);
        assert(trivial_pow(2, 10) == 1024);
        assert(trivial_pow(2, 20) == 1048576);
        release_memory(&local_buffer);

        assert(pow_by_squaring(0, 0) == 1);
        assert(pow_by_squaring(0, 1) == 0);
        assert(pow_by_squaring(0, 54) == 0);
        assert(pow_by_squaring(0, 5554135) == 0);
        assert(pow_by_squaring(2, 0) == 1);
        assert(pow_by_squaring(1048576, 0) == 1);
        assert(pow_by_squaring(5465313, 1) == 5465313);
        assert(pow_by_squaring(5465313, 2) == 5465313ull*5465313ull);
        assert(pow_by_squaring(2, 10) == 1024);
        assert(pow_by_squaring(2, 20) == 1048576);
        assert(pow_by_squaring(3, 10) == 59049);
        assert(pow_by_squaring(3, 20) == 59049ull*59049ull);
        assert(pow_by_squaring(9, 10) == 59049ull*59049ull);
        assert(pow_by_squaring(5465, 3) == 5465ull*5465ull*5465ull);
    }
    
    runtime_proc run_typed_tests(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        test_misc(generator, random_runs);
    }

    template <typename T>
    runtime_proc run_untyped_tests(Memory_Resource* upstream, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        test_add_overflow<T>(upstream, generator, random_runs);
        test_sub_overflow<T>(upstream, generator, random_runs);
        test_complement_overflow<T>(upstream, generator, random_runs);
        test_shift_overflow<T>(upstream, generator, random_runs);
        test_mul_overflow<T>(upstream, generator, random_runs, controlled_runs);
        test_div_overflow_low<T>(upstream, generator, random_runs, controlled_runs);
        test_mul_quadratic<T>(upstream, generator, random_runs);
        test_div_bit_by_bit<T>(upstream, generator, random_runs);
        test_pow_by_squaring<T>(upstream, generator, random_runs, controlled_runs);

        
        const size_t quarter_runs = random_runs / 4;
        if constexpr(sizeof(T) >= sizeof(u8))
            test_to_base<T, u8>(upstream, generator, quarter_runs);

        if constexpr(sizeof(T) >= sizeof(u16))
            test_to_base<T, u16>(upstream, generator, quarter_runs);

        if constexpr(sizeof(T) >= sizeof(u32))
            test_to_base<T, u32>(upstream, generator, quarter_runs);

        if constexpr(sizeof(T) >= sizeof(u64))
            test_to_base<T, u64>(upstream, generator, quarter_runs);
        
    }

    runtime_proc run_tests()
    {
        std::random_device os_seed;
        const u32 seed = os_seed();

        Random_Generator generator(seed);
        const size_t random_runs = 10000;
        const size_t controlled_runs = 500;

        //Without: 3592801800
        //With:    3108090100
        for(int i = 0; i < 10; i++)
        {
            std::cout << "ellapsed: " << ellapsed_time([&](){
                //std::pmr::unsynchronized_pool_resource pool;
                //Memory_Resource* res = &pool;

                #ifdef USE_CUSTOM_LIB
                //jot::Unbound_Arena_Resource arena;
                Memory_Resource* res = jot::new_delete_resource();
                #else
                Memory_Resource* res = std::pmr::new_delete_resource();
                #endif
                run_typed_tests(res, &generator, random_runs, controlled_runs);
                run_untyped_tests<u8>(res, &generator, random_runs, controlled_runs);
                run_untyped_tests<u16>(res, &generator, random_runs, controlled_runs);
                run_untyped_tests<u32>(res, &generator, random_runs, controlled_runs);
                run_untyped_tests<u64>(res, &generator, random_runs, controlled_runs);
            }) << " ms" << std::endl;
        }
    }
}