#pragma once

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <cctype>
#include <type_traits>

#include <iostream>
#include <random>
#include <limits>
#include <initializer_list>
#include <memory_resource>
#include <vector>
#include <array>
#include <algorithm>

#define USE_CUSTOM_LIB

#include "jot/meta.h"

#ifdef USE_CUSTOM_LIB
#include "jot/stack.h"
#include "jot/allocator_stack.h"
#include "jot/allocator_ring.h"
#include "jot/defer.h"
#endif // USE_CUSTOM_LIB

#include "benchmark.h"
#include "tiny_big_num.h"

#pragma once


#define let const auto
#define mut auto
#define proc constexpr auto
#define func [[nodiscard]] constexpr auto
#define runtime_proc [[nodiscard]] auto
#define runtime_func auto

#define cast(...) (__VA_ARGS__)

namespace tiny_num::test
{
    using u8 = std::uint8_t;
    using u16 = std::uint16_t;
    using u32 = std::uint32_t;
    using u64 = std::uint64_t;

    using i8 = std::int8_t;
    using i16 = std::int16_t;
    using i32 = std::int32_t;
    using i64 = std::int64_t;

    using b8 = bool;
    using b16 = std::uint16_t;
    using b32 = std::uint32_t;
    using b64 = std::uint64_t;

    using f32 = float;
    using f64 = double;

    using cstring = const char*;
    using isize = i64;
    using usize = size_t;

    using umax = u64;

    using std::move;
    using std::size;
    using std::data;

    #ifdef USE_CUSTOM_LIB
        template <typename T>
        using Vector = jot::Stack<T>;

        using Memory_Resource = jot::Allocator;

        template<typename T>
        func make_sized_vector(size_t size, Memory_Resource* resource) -> Vector<T>
        {
            using namespace jot;
            Vector<T> vec = {resource};
            if(size != 0)
                force(resize(&vec, cast(isize) size));
            return vec;
        }

        template<typename T>
        func push_back(Vector<T>* vec, T pushed)
        {
            jot::force(jot::push(vec, move(pushed)));
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

        template<typename T>
        func make_sized_vector(size_t size, Memory_Resource* resource) -> Vector<T>
        {
            Vector<T> vec(std::move(resource));
            if(size != 0)
                vec.resize(size);
            return vec;
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
        copy_slice<T>(&vec_slice, slice, Iter_Direction::NO_ALIAS);

        return vec;
    }

    template<typename T>
    func make_vector_of_digits(umax val, Memory_Resource* resource) -> Vector<T>
    {
        constexpr size_t count = digits_to_represent<T, umax>();
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
        size_t prefix_size = 0;
        size_t content_size = 0;
        size_t postfix_size = 0;
    };

    template <typename T>
    func make_padded_vector_of_digits(umax val, Memory_Resource* resource, size_t pref_size = 1, size_t post_size = 1) -> Padded_Vector<T>
    {
        constexpr size_t count = digits_to_represent<T, umax>();
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
        umax value;
        umax overflow;

        bool constexpr operator ==(Batch_Op_Result const&) const noexcept = default;
    };

    template <typename T>
    constexpr size_t MAX_TYPE_SIZE_FRACTION = sizeof(umax) / sizeof(T);

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

    template <typename T, typename ReturnT = umax>
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

    u64 single_pow_by_squaring_native(u64 x, u64 n)
    {
        if (x <= 1) 
            return x;

        u64 y = 1;
        for (; n != 0; n >>= 1, x *= x)
            if (n & 1)
                y *= x;

        return y;
    }

    u64 single_pow_trivial(u64 x, u64 n)
    {
        u64 res = 1;
        while(n > 0)
        {
            res *= x;
            n--;
        }

        return res;
    }

    u64 single_log_heyley(u64 of_value, u64 base)
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
            u64 value_from_curr_approximate = single_pow_by_squaring(b, y);
            u64 c_x = value_from_curr_approximate;

            prev_y = y;
            if(c_x <= x)
                break;

            u64 new_y = y - 2*(c_x - x) / (x + c_x);

            //if we didnt move at all we are one above the answer
            if(new_y == y)
                return new_y - 1;

            y = new_y;
        }

        return prev_y;

    }

    u64 single_log_bsearch(u64 of_value, u64 base)
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

        while(true)
        {
            u64 mid = (curr_lo + curr_hi) / 2;
            u64 curr_approx = single_pow_by_squaring(b, mid);

            if(curr_approx > x)
                curr_hi = mid;
            else if(curr_approx < x)
                curr_lo = mid + 1;
            else
                return mid;

            if(curr_hi <= curr_lo)
                break;
        }

        return curr_lo - 1;
    }

    runtime_proc test_misc(Random_Generator* generator, size_t random_runs)
    {
        using umax = umax;

        assert(find_last_set_bit<umax>(0) == 0);
        assert(find_last_set_bit<umax>(1) == 0);
        assert(find_last_set_bit<umax>(2) == 1);
        assert(find_last_set_bit<umax>(3) == 1);
        assert(find_last_set_bit<umax>(4) == 2);
        assert(find_last_set_bit<umax>(7) == 2);
        assert(find_last_set_bit<umax>(8) == 3);
        assert(find_last_set_bit<umax>(76) == 6);
        assert(find_last_set_bit<umax>(513) == 9);
        assert(find_last_set_bit<umax>(641) == 9);
        assert(find_last_set_bit<umax>(1024) == 10);
        assert(find_last_set_bit<u64>(0xf4f1f63db9b7e800) == 63);
        assert(find_last_set_bit<umax>(0b010100) == 4);
        assert(find_last_set_bit<umax>(0b010000) == 4);

        assert(find_first_set_bit<umax>(0) == 0);
        assert(find_first_set_bit<umax>(0b000001) == 0);
        assert(find_first_set_bit<umax>(0b010100) == 2);
        assert(find_first_set_bit<umax>(0b010101) == 0);
        assert(find_first_set_bit<umax>(0b010000) == 4);
        assert(find_first_set_bit<umax>(0b010100) == 2);
        assert(find_first_set_bit<umax>(0b010010) == 1);

        assert(pop_count<umax>(0b0000'0000) == 0);
        assert(pop_count<umax>(0b0000'0001) == 1);
        assert(pop_count<umax>(0b1000'0000) == 1);
        assert(pop_count<umax>(0b1010'0000) == 2);
        assert(pop_count<umax>(0b1010'0011) == 4);
        assert(pop_count<umax>(0b0111'1111) == 7);
        assert(pop_count<umax>(cast(u32) -1) == 32);
        assert(pop_count<umax>(cast(umax) -1) == 64);

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


        assert(single_pow_by_squaring(3,3) == 27);

        assert(single_root_shifting(8, 3) == 2);
        assert(single_root_shifting(9, 3) == 2);
        assert(single_root_shifting(27, 3) == 3);
        assert(single_root_shifting(28, 3) == 3);
        assert(single_root_shifting(29, 3) == 3);
        assert(single_root_shifting(34, 5) == 2);
        assert(single_root_shifting(15625, 3) == 25);
        assert(single_root_shifting(0, 0) == 1);
        assert(single_root_shifting(1, 0) == 1);
        assert(single_root_shifting(4096, 6) == 4);
        assert(single_root_shifting(18446744073709551614, 17) == 13);

        assert(single_root_newton(8, 3) == 2);
        assert(single_root_newton(9, 3) == 2);
        assert(single_root_newton(27, 3) == 3);
        assert(single_root_newton(28, 3) == 3);
        assert(single_root_newton(29, 3) == 3);
        assert(single_root_newton(34, 5) == 2);
        assert(single_root_newton(15625, 3) == 25);
        assert(single_root_newton(0, 0) == 1);
        assert(single_root_newton(1, 0) == 1);
        assert(single_root_newton(4096, 6) == 4);


        assert(single_log_heyley(25, 2) == 4);
        assert(single_log_heyley(255, 2) == 7);
        assert(single_log_heyley(256, 2) == 8);
        assert(single_log_heyley(4, 4) == 1);
        assert(single_log_heyley(9, 3) == 2);
        assert(single_log_heyley(26, 3) == 2);
        assert(single_log_heyley(27, 3) == 3);
        assert(single_log_heyley(64, 3) == 3);
        assert(single_log_heyley(81, 3) == 4);
        assert(single_log_heyley(16, 4) == 2);
        assert(single_log_heyley(17, 4) == 2);
        assert(single_log_heyley(27, 4) == 2);
        assert(single_log_heyley(63, 4) == 2);
        assert(single_log_heyley(64, 4) == 3);


        assert(single_log_bsearch(25, 2) == 4);
        assert(single_log_bsearch(255, 2) == 7);
        assert(single_log_bsearch(256, 2) == 8);
        assert(single_log_bsearch(4, 4) == 1);
        assert(single_log_bsearch(9, 3) == 2);
        assert(single_log_bsearch(26, 3) == 2);
        assert(single_log_bsearch(27, 3) == 3);
        assert(single_log_bsearch(64, 3) == 3);
        assert(single_log_bsearch(81, 3) == 4);
        assert(single_log_bsearch(16, 4) == 2);
        assert(single_log_bsearch(17, 4) == 2);
        assert(single_log_bsearch(27, 4) == 2);
        assert(single_log_bsearch(63, 4) == 2);
        assert(single_log_bsearch(64, 4) == 3);

        assert((highest_power_of_in_type(2, 2) == Power_And_Powed{1, 2}));
        assert((highest_power_of_in_type(1, 2) == Power_And_Powed{0, 1}));
        assert((highest_power_of_in_type(16, 2) == Power_And_Powed{15, 1 << 15}));
        assert((highest_power_of_in_type(8, 3) == Power_And_Powed{5, 243}));

        assert((highest_power_of_in_type(0, 10) == Power_And_Powed{0, 1}));
        assert((highest_power_of_in_type(1, 10) == Power_And_Powed{0, 1}));
        assert((highest_power_of_in_type(2, 10) == Power_And_Powed{0, 1}));
        assert((highest_power_of_in_type(4, 10) == Power_And_Powed{1, 10}));
        assert((highest_power_of_in_type(8, 10) == Power_And_Powed{2, 100}));
        assert((highest_power_of_in_type(16, 10) == Power_And_Powed{4, 10'000}));
        assert((highest_power_of_in_type(32, 10) == Power_And_Powed{9, 1'000'000'000}));
        assert((highest_power_of_in_type(33, 10) == Power_And_Powed{9, 1'000'000'000}));
        assert((highest_power_of_in_type(34, 10) == Power_And_Powed{10, 10'000'000'000}));
        assert((highest_power_of_in_type(64, 10) == Power_And_Powed{19, 10'000'000'000'000'000'000}));

        constexpr u64 max = std::numeric_limits<u64>::max()/2;
        constexpr u64 min = 2;
        let val_dist = make_exponential_distribution<u64>(min, max);

        assert(cast(u64) 1 << 42 == 4398046511104);
        assert(single_log_heyley(4398046511104, 2) == 42);

        for(size_t i = 0; i < random_runs; i++)
        {
            u64 val = cast(u64) val_dist(*generator);
            let power_dist = make_exponential_distribution<u64>(min, val);
            u64 power = cast(u64) power_dist(*generator);

            u64 rooted = single_root_shifting(val, power);
            u64 powed_low = single_pow_by_squaring(rooted, power);
            u64 back_down = single_root_shifting(powed_low, power);
            if(rooted > 1 && powed_low < 43980465111) //so that no overflow
            {
                u64 loged_heyley = single_log_heyley(powed_low, rooted);
                u64 loged_bsearch = single_log_heyley(powed_low, rooted);
                assert(loged_heyley == power);
                assert(loged_bsearch == power);
            }

            if(power < 100)
            {
                u64 ref_powed_low = single_pow_trivial(rooted, power);
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


            u64 rooted_again = single_root_shifting(rooted, power);
            u64 rooted_again_newton = single_root_newton(rooted, power);

            assert(rooted_again == rooted_again_newton);
        }
    }

    template <typename T>
    runtime_proc test_add_overflow(Memory_Resource* memory, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        runtime_proc test_add_overflow_batch = [&](umax left_, umax right_, umax carry_in, Res expected, bool check_overflow = true) -> bool {
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Padded pad_right = make_padded_vector_of_digits<T>(right_, resource, 1);

            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);
            Slice pad_right_s = to_slice<T>(&pad_right.vector);

            let loc = Op_Location::OUT_OF_PLACE;

            let res1 = ::batch_add_overflow_long<T>(&out_s,  to_slice(left), to_slice(right), loc, cast(T) carry_in);
            let res2 = ::batch_add_overflow_long<T>(&pad_left_s, content(pad_left), to_slice(right), loc, cast(T) carry_in); //test self assign
            let res3 = ::batch_add_overflow_long<T>(&pad_right_s, to_slice(left), content(pad_right), loc, cast(T) carry_in); 

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
    
        runtime_proc auto_test_add = [&](umax left, umax right, umax carry_in = 0) -> bool {
            umax highest_one = cast(umax) 1 << (BIT_SIZE<umax> - 1);
            umax adjusted_left = left & ~highest_one;
            umax adjusted_right = right & ~highest_one;
            umax adjusted_carry = cast(umax) low_bits(cast(T) carry_in);
            umax expected = single_add_no_overflow(adjusted_left, adjusted_right, adjusted_carry);

            return test_add_overflow_batch(adjusted_left, adjusted_right, adjusted_carry, Res{expected, 0}, false);
        };

        runtime_proc test_single_add = [&](umax left, umax single){
            if(left == 0 || single == 0)
                return true;

            umax highest_one = cast(umax) 1 << (BIT_SIZE<umax> - 1);
            umax adjusted_left = left & ~highest_one;
            T adjusted_single = cast(T) (single & ~highest_one);

            umax expected = adjusted_left + adjusted_single;

            mut left1 = make_vector_of_digits<T>(adjusted_left, resource); 
            mut left2 = make_vector_of_digits<T>(adjusted_left, resource);

            Slice left1_s = to_slice(&left1);
            Slice left2_s = to_slice(&left2);

            let res2 = batch_add_overflow_long<T>(&left1_s, left1_s, adjusted_single, Op_Location::IN_PLACE);
            let res1 = batch_add_overflow_long<T>(&left2_s, left2_s, Slice{&adjusted_single, 1}, Op_Location::IN_PLACE);

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
            umax left = distribution(*generator);
            umax right = distribution(*generator);
            umax carry = distribution(*generator);
            assert(auto_test_add(left, right, carry));
            assert(test_single_add(left, cast(T) right));
        };
    }

    
    template <typename T>
    runtime_proc test_sub_overflow(Memory_Resource* memory, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        runtime_proc test_sub_overflow_batch = [&](umax left_, umax right_, umax carry_in, Res expected, bool check_overflow = true) -> bool {
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Padded pad_right = make_padded_vector_of_digits<T>(right_, resource, 1);

            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);

            let loc = Op_Location::OUT_OF_PLACE;

            let res1 = ::batch_sub_overflow_long<T>(&out_s,  to_slice(left), to_slice(right), loc, cast(T) carry_in);
            let res2 = ::batch_sub_overflow_long<T>(&pad_left_s, content(pad_left), to_slice(right), loc, cast(T) carry_in); //test self assign

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

        runtime_proc auto_test_sub = [&](umax left, umax right, umax carry_in = 0) -> bool {
            if(left < right)
                std::swap(left, right);

            umax highest_one = cast(umax) 1 << (BIT_SIZE<umax> - 1);
            umax adjusted_right = right & ~highest_one;
            umax adjusted_carry = carry_in > low_mask<umax>();

            umax adjusted_left = left;
            if(left < adjusted_carry + adjusted_right)
                adjusted_left = single_add_no_overflow(left, adjusted_carry);

            umax expected = adjusted_left - adjusted_right - adjusted_carry;
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
            umax left = distribution(*generator);
            umax right = distribution(*generator);
            umax carry = distribution(*generator);
            assert(auto_test_sub(left, right, carry));
        };
    }

    template <typename T>
    runtime_proc test_complement_overflow(Memory_Resource* memory, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        //we dont test oveflow of this op 
        runtime_proc batch_complement_overflow = [&](umax left_, umax carry_in = 1) -> umax{
            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Slice left_s = content(&pad_left);
            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);

            let loc = Op_Location::OUT_OF_PLACE;

            let res1 = ::batch_complement_overflow<T>(&out_s, left_s, cast(T) carry_in);
            let res2 = ::batch_complement_overflow<T>(&pad_left_s, left_s, cast(T) carry_in); //test self assign

            let num1 = unwrap(to_number(res1.slice));
            let num2 = unwrap(to_number(res2.slice));

            assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
            assert(num1 == num2 && "self assign must work correctly");

            return num1;
        };


        runtime_proc auto_test_complement = [&](umax left, umax carry_in = 1) -> bool{
            umax complemented = ~left + carry_in;
            umax cropped_complemented = complemented;

            size_t bits = 0;
            if(left == 0)
                bits = 0;
            else if(BIT_SIZE<T> >= BIT_SIZE<umax>)
                bits = BIT_SIZE<umax>;
            else
            {
                for(umax temp_left = left; temp_left > 0;)
                {
                    temp_left >>= BIT_SIZE<T>;
                    bits += BIT_SIZE<T>;
                }
            }
        
            if(bits < BIT_SIZE<umax>)
                cropped_complemented = low_bits(complemented, bits);

            return batch_complement_overflow(left, carry_in) == cropped_complemented;
        };
    
        if constexpr(std::is_same_v<T, u8>)
        {
            assert(batch_complement_overflow(0) == 0);
            assert(batch_complement_overflow(1) == cast(u8)(-1));
            assert(batch_complement_overflow(5) == cast(u8)(-5));
            assert(batch_complement_overflow(0xFF) == cast(u8)(-0xFF));
            assert(batch_complement_overflow(0xABCD) == cast(u16)(-0xABCD));
            assert(batch_complement_overflow(0xFFFF0000) == cast(u32)(-cast(i64)0xFFFF0000));
        }

        let exponential = make_exponential_distribution<T, umax>();
        let uniform = std::uniform_int_distribution<>(0, 1);
        for(size_t i = 0; i < random_runs; i++)
        {
            umax left = exponential(*generator);
            umax carry = uniform(*generator);
            assert(auto_test_complement(left, carry));
        };
    }

    template <typename T>
    runtime_proc test_shift_overflow(Memory_Resource* memory, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        runtime_proc shift_both = [&](umax left_, umax right, umax carry_in, bool up_down, Iter_Direction direction) -> Batch_Op_Result{
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
                res1 = ::batch_shift_up_overflow<T>(&out_s, left_s, right, direction, cast(T) carry_in);
                res2 = ::batch_shift_up_overflow<T>(&pad_left_s, left_s, right, direction, cast(T) carry_in);
            }
            else
            {
                res1 = ::batch_shift_down_overflow<T>(&out_s, left_s, right, direction, cast(T) carry_in);
                res2 = ::batch_shift_down_overflow<T>(&pad_left_s, left_s, right, direction, cast(T) carry_in);
            }

            let num1 = unwrap(to_number(res1.slice));
            let num2 = unwrap(to_number(res2.slice));

            assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
            assert(res1.overflow == res2.overflow);
            assert(num1 == num2 && "self assign must work correctly");

            return Batch_Op_Result{num1, res1.overflow};
        };

        runtime_proc batch_shift_up_overflow = [&](umax left, umax right, umax carry_in = 0) -> Batch_Op_Result{
            let forward = shift_both(left, right, carry_in, true, Iter_Direction::FORWARD);
            let backward = shift_both(left, right, carry_in, true, Iter_Direction::BACKWARD);

            assert(forward.value == backward.value);
            assert(forward.overflow == backward.overflow);

            return forward;
        };

        runtime_proc batch_shift_down_overflow = [&](umax left, umax right, umax carry_in = 0) -> Batch_Op_Result{
            let forward = shift_both(left, right, carry_in, false, Iter_Direction::FORWARD);
            let backward = shift_both(left, right, carry_in, false, Iter_Direction::BACKWARD);

            assert(forward.value == backward.value);
            assert(forward.overflow == backward.overflow);

            return forward;
        };


        if constexpr(std::is_same_v<T, u8>)
        {
            assert((batch_shift_up_overflow(0, 0, 0) == Res{0, 0}));
            assert((batch_shift_up_overflow(0, 1, 0) == Res{0, 0}));
            assert((batch_shift_up_overflow(0, 1, 0xFF) == Res{0, 0xFF})); //should be consistent with multiply

            assert((batch_shift_up_overflow(1, 1) == Res{2, 0}));
            assert((batch_shift_up_overflow(1, 4) == Res{16, 0}));
            assert((batch_shift_up_overflow(1, 4, 0b1010'0000) == Res{0b11010, 0}));
            assert((batch_shift_up_overflow(11, 4) == Res{176, 0}));
            assert((batch_shift_up_overflow(0xAB, 0x4) == Res{0xB0, 0xA}));
            assert((batch_shift_up_overflow(0xAB, 0x4, 0xD0) == Res{0xBD, 0xA}));
            assert((batch_shift_up_overflow(0xAB, 0x4, 0x0D) == Res{0xB0, 0xA}));
            assert((batch_shift_up_overflow(0b0101'0101'0111'0001, 7) == Res{0b1011'1000'1000'0000, 0b010'1010}));
            assert((batch_shift_up_overflow(0b0101'0101'0111'0001, 0) == Res{0b0101'0101'0111'0001, 0}));
            assert((batch_shift_up_overflow(0b0101'0101'0111'0001, 7, 0b11001100) == Res{0b1011'1000'1110'0110, 0b010'1010}));

            assert((batch_shift_down_overflow(0, 0, 0) == Res{0, 0}));
            assert((batch_shift_down_overflow(0, 1, 0) == Res{0, 0}));
            assert((batch_shift_down_overflow(0, 1, 0xFF) == Res{0, 0xFF}));
            assert((batch_shift_down_overflow(1, 1) == Res{0, 0x80}));
            assert((batch_shift_down_overflow(0b10, 1) == Res{0b1, 0}));
            assert((batch_shift_down_overflow(0b110000, 3) == Res{0b110, 0}));
            assert((batch_shift_down_overflow(0b0101'0011, 1) == Res{0b101001, 0x80}));
            assert((batch_shift_down_overflow(0b0101'0011, 4) == Res{0b101, 0b0011'0000}));
            assert((batch_shift_down_overflow(0b0101'0011, 4, 0b0000'1111) == Res{0b1111'0101, 0b0011'0000}));
            assert((batch_shift_down_overflow(0b0101'0011, 0) == Res{0b0101'0011, 0}));
            assert((batch_shift_down_overflow(0b0101'0011, 7) == Res{0b0, 0b10100110}));
            assert((batch_shift_down_overflow(0b0101'0011, 7, 0b1100'0111) == Res{0b1000'1110, 0b10100110}));
            assert((batch_shift_down_overflow(0xFFEEDD, 4) == Res{0x0FFEED, 0xD0}));
            assert((batch_shift_down_overflow(0xFFEEDD, 4, 0x0A) == Res{0xAFFEED, 0xD0}));
            assert((batch_shift_down_overflow(0b0101010101110001, 5) == Res{0b0000001010101011, 0b10001000}));
            assert((batch_shift_down_overflow(0b0101010101110001, 5, 0b1010'1111) == Res{0b0111101010101011, 0b10001000}));
        }
    }

    struct Adjusted
    {
        umax left;
        umax right;
    };

    runtime_func adjust_to_not_mul_overflow(umax left, umax right) -> Adjusted
    {
        using umax = umax;
        // (by doing find_last_set_bit) and then shifting it so that they will not overflow

        umax adjusted_left = left;
        umax adjusted_right = right;

        size_t log_left = find_last_set_bit(left) + 1; 
        size_t log_right = find_last_set_bit(right) + 1; 
        //find_last_set_bit returns position of highest set bit but we want "under which power of two is the whole number" 
        // => we add 1

        size_t combined_log = log_left + log_right;

        //if can overflow we shift down so that its not possible
        if(combined_log > BIT_SIZE<umax>)
        {
            size_t diff = combined_log - BIT_SIZE<umax>;
            size_t shift = (diff + 1) / 2;

            adjusted_left >>= shift;
            adjusted_right >>= shift;
        }

        assert(find_last_set_bit(adjusted_left) + find_last_set_bit(adjusted_right) + 2 <= BIT_SIZE<umax>);
        return {adjusted_left, adjusted_right};
    }

    template <typename T>
    runtime_proc test_mul_overflow(Memory_Resource* memory, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        runtime_proc test_mul_overflow_batch = [&](umax left_, umax right, umax carry_in, Res expected, Optim_Info optims, bool check_overflow = true) -> bool {
            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Slice left_s = content(&pad_left);
            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = to_slice<T>(&pad_left.vector);

            let res1 = ::batch_mul_overflow<T>(&out_s, left_s, cast(T) right, optims, cast(T) carry_in);
            let res2 = ::batch_mul_overflow<T>(&pad_left_s, left_s, cast(T) right, optims, cast(T) carry_in);

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

        runtime_proc auto_test_mul = [&](umax left, umax right, Optim_Info optims) -> bool
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

                assert((single_mul_overflow<T>(0xFF, 0x10, 0xC) == Overflow<T>{0xFC, 0xF}));
                assert(test_mul_overflow_batch(0xFFEE, 0x80, 0, Res{0xF700, 0x7F}, optims));
                assert(test_mul_overflow_batch(0xFFEE, 0x80, 0, Res{0xF700, 0x7F}, optims));
                assert(test_mul_overflow_batch(0xA1, 0x11, 0, Res{0xB1, 0xA}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 1, 0, Res{0xABCDEF00, 0}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 0, 0, Res{0, 0}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 0xAB, 0, Res{0xC28EA500, 0x72}, optims));
                assert(test_mul_overflow_batch(13745, 137, 0, Res{0xBBB9, 0x1C}, optims));
            }
        }

        assert(auto_test_mul(0x14156f, 0x3c, Optim_Info{}));

        let exponential = make_exponential_distribution<T, umax>();
        let uniform = std::uniform_int_distribution<long long>(0, FULL_MASK<T> >> 1);
        for(size_t i = 0; i < random_runs; i++)
        {
            umax left = exponential(*generator);
            umax right = uniform(*generator);
            let optims = generate_random_optims(generator);
            assert(auto_test_mul(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_div_overflow_low(Memory_Resource* memory, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        assert((single_div_overflow<T>(0x12, 0x0A, 0) == Overflow<T>{1, 0x08}));

        Memory_Resource* resource = memory;

        runtime_proc test_div_overflow_low_batch = [&](umax left_, umax right, umax carry_in, Res expected, Optim_Info optims) -> bool {
            constexpr size_t padding_before = 2;

            Padded pad_left = make_padded_vector_of_digits<T>(left_, resource, 1);
            Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1, resource);

            Slice left_s = content(&pad_left);
            Slice out_s = to_slice<T>(&out);
            Slice pad_left_s = slice(to_slice<T>(&pad_left.vector), 1);

            let res1 = ::batch_div_overflow<T>(&out_s, left_s, cast(T) right, optims, cast(T) carry_in);
            let res2 = ::batch_div_overflow<T>(&pad_left_s, left_s, cast(T) right, optims, cast(T) carry_in);

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
    
        runtime_proc auto_test_div = [&](umax left, umax right, Optim_Info optims) -> bool {
            if(right == 0)
                return true;
        
            umax adjusted_right = low_bits(cast(T) right);
            umax res = left / adjusted_right;
            umax rem = left % adjusted_right;

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
            }
        }

        let exponential = make_exponential_distribution<T, umax>();
        let uniform = std::uniform_int_distribution<long long>(0, low_mask<T>());
        for(size_t i = 0; i < random_runs; i++)
        {
            umax left = exponential(*generator);
            umax right = uniform(*generator);
            let optims = generate_random_optims(generator);
            assert(auto_test_div(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_mul_quadratic(Memory_Resource* memory, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        runtime_proc test_fused_mul_add = [&](umax left_, umax right_, umax coef, umax expected, Optim_Info optims = Optim_Info{}) -> bool{
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);
            Vector output = make_sized_vector<T>(size(left), resource);

            Slice left_s = to_slice(&left);
            Slice right_s = to_slice(&right);
            Slice output_s = to_slice(&output);

            let res_out = ::batch_fused_mul_add_overflow<T>(&output_s, left_s, right_s, cast(T) coef, optims, Op_Location::OUT_OF_PLACE);
            let res_in = ::batch_fused_mul_add_overflow<T>(&left_s, left_s, right_s, cast(T) coef, optims, Op_Location::IN_PLACE);
        
            push_back(&output, res_out.overflow);
            push_back(&left, res_in.overflow);

            left_s = to_slice(&left);
            output_s = to_slice(&output);

            let num_out = unwrap(to_number(left_s));
            let num_in = unwrap(to_number(output_s));

            return expected == num_out && expected == num_in; 
        };

        runtime_proc auto_test_fused_mul_add = [&](umax left, umax right, umax coef, Optim_Info optims = Optim_Info{}) -> bool{
            return test_fused_mul_add(left, right, coef, left + coef*right, optims);
        };

        runtime_proc test_mul_quadratic = [&](umax left_, umax right_, umax expected, Optim_Info optims = Optim_Info{}) -> bool{
            Vector left = make_vector_of_digits<T>(left_, resource);
            Vector right = make_vector_of_digits<T>(right_, resource);

            size_t to_size = required_mul_out_size(size(left), size(right));
            size_t aux_size1 = required_mul_quadratic_aux_size(size(left), size(right));
            size_t aux_size2 = required_mul_karatsuba_aux_size(size(left), size(right));

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

        runtime_proc auto_test_mul_quadratic = [&](umax left, umax right, Optim_Info optims) -> bool
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

        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x38, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0xFFFF, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x1010, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x7838, shift));
        assert(auto_test_mul_quadratic(0x1d15a574ce80, 0x27838, shift));
        assert(auto_test_mul_quadratic(60136640465, 262154, shift));

        let distribution = make_exponential_distribution<T>();
        for(size_t i = 0; i < random_runs; i++)
        {
            Optim_Info optims = generate_random_optims(generator);
            umax left = distribution(*generator);
            umax right = distribution(*generator);

            assert(auto_test_mul_quadratic(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_div_bit_by_bit(Memory_Resource* memory, Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        enum Will_Fail
        {
            FAIL,
            OKAY
        };

        Memory_Resource* resource = memory;

        runtime_proc test_div_bit_by_bit = [&](umax num_, umax den_, umax ex_quotient, umax ex_remainder, Will_Fail expected_fail = OKAY, Optim_Info optims = Optim_Info{}) -> bool{
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


        runtime_proc auto_test_div_bit_by_bit = [&](umax left, umax right, Optim_Info optims = Optim_Info{}) -> bool{
            if(right == 0)
                return test_div_bit_by_bit(left, right, 0, 0, FAIL);

            umax quotient = left / right;
            umax remainder = left % right;
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
            Optim_Info optims = generate_random_optims(generator);
            umax left = distribution(*generator);
            umax right = distribution(*generator);
            assert(auto_test_div_bit_by_bit(left, right, optims));
        };
    }

    template <typename T>
    void println(Slice<T> slice)
    {
        std::cout << "[";
        if(slice.size > 0)
        {
            std::cout << cast(umax) slice[0];

            for(size_t i = 1; i < slice.size; i++)
                std::cout << ", " << cast(umax) slice[i];
        }
        std::cout << "]\n";
    }

    template <typename Num, typename Rep, typename Conv_Fn>
    func to_base(Slice<const Num> num, Num base, Conv_Fn conv, Optim_Info optims, Memory_Resource* resource) -> Vector<Rep> {
        Vector<Num> temp = make_vector_of_slice<Num>(num, resource);
        Vector<Rep> out = make_sized_vector<Rep>(required_to_base_out_size<Num>(num.size, base), resource);

        Slice<Num> temp_s = to_slice(&temp);
        Slice<Rep> out_s = to_slice(&out);

        Slice<Rep> converted = (to_base<Num, Rep>(&out_s, &temp_s, num, base, conv, optims));
        assert(is_striped_representation(converted));

        resize(&out, converted.size);
        return out;
    };

    template <typename Num, typename Rep, typename Conv_Fn>
    func from_base(Slice<const Rep> rep, Num base, Conv_Fn conv, Optim_Info optims, Memory_Resource* resource) -> Vector<Num> {
        size_t required_size = required_from_base_out_size<Num>(rep.size, base);
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

    template <typename Num = umax>
    func to_base(Slice<const Num> num, Num base, Memory_Resource* resource) -> Vector<char> {
        assert(base > 2);
        assert(base <= 36);
        let conversion = [&](Num value) -> char {

            constexpr char val_to_char_table[36] = {
                '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
                'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                'U', 'V', 'W', 'X', 'Y', 'Z'
            };

            assert(value < 36);
            return val_to_char_table[value];
        };

        return to_base<Num, char>(num, base, conversion, Optim_Info{}, resource);
    }

    template <typename Num = umax>
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
    func native_to_base(umax num, Num base, Memory_Resource* resource) -> Vector<Rep> 
    {
        using umax = umax;
        constexpr size_t num_digits = digits_to_represent<Num, umax>();
        const size_t max_size = required_to_base_out_size<Num>(num_digits, base);

        Vector<Rep> converted = make_sized_vector<Rep>(max_size, resource);
        size_t size = 0;
        for(umax curr_val = num; curr_val != 0; curr_val /= base, size ++)
        {
            umax rem = curr_val % base;
            converted[size] = cast(Rep) rem;
        }

        resize(&converted, size);
        Slice<Rep> converted_s = to_slice(&converted);
        reverse(&converted_s);

        return converted;
    };

    template <typename Num, typename Rep>
    func native_from_base(Slice<const Rep> rep, Num base) -> umax 
    {
        using umax = umax;
        assert(is_striped_representation(rep));

        umax value = 0;
        for(size_t i = 0; i < rep.size; i++)
        {
            umax muled = value * base;
            umax added = muled + rep[i];
            assert(muled >= value);
            assert(added >= muled);

            value = added;
        }

        return value;
    };

    template <typename Num, typename Rep>
    runtime_proc test_to_base(Memory_Resource* memory, Random_Generator* generator, size_t random_runs)
    {
        static_assert(sizeof(Num) >= sizeof(Rep));
        using umax = umax;

        Memory_Resource* resource = memory;

        proc id_to_conversion = [](Num num) -> Rep {return cast(Rep) num;};
        proc id_from_conversion = [](Rep num) -> Num {return cast(Num) num;};

        runtime_proc test_to_base = [&](Slice<const Num> num, Num base, Slice<const Rep> expected_rep, Optim_Info optims) -> bool 
        {
            Vector<Rep> converted = to_base<Num, Rep>(num, base, id_to_conversion, optims, resource);
            assert(is_striped_representation(expected_rep));

            return is_equal<Rep>(to_slice(converted), expected_rep);
        };

        runtime_proc manual_test_to_base = [&](umax num, Num base, std::initializer_list<Rep> expected_rep, Optim_Info optims = Optim_Info{}) -> bool 
        {
            Vector<Num> num_ = make_vector_of_digits<Num>(num, resource);
            Slice<const Rep> expected_rep_s = {std::data(expected_rep), std::size(expected_rep)};
            return test_to_base(to_slice(num_), base, expected_rep_s, optims);
        };

        runtime_proc auto_test_to_base = [&](umax num, Num base, Optim_Info optims) -> bool 
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

        runtime_proc manual_test_from_base = [&](std::initializer_list<Rep> rep, Num base, umax expected_num, Optim_Info optims = Optim_Info{}) -> bool 
        {
            Vector<Num> expected_num_ = make_vector_of_digits<Num>(expected_num, resource);
            Slice<const Rep> rep_s = {std::data(rep), std::size(rep)};
            return test_from_base(rep_s, base, to_slice(expected_num_), optims);
        };

        runtime_proc auto_test_from_base = [&](Slice<const Rep> rep, Num base, Optim_Info optims) -> bool 
        {
            umax num = native_from_base<Num, Rep>(rep, base);
            Vector<Num> expected_num = make_vector_of_digits<Num>(num, resource);
            return test_from_base(rep, base, to_slice(expected_num), optims);
        };

        runtime_proc test_to_and_fro = [&](umax value, Num base, Optim_Info optims) -> bool
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
        assert(auto_test_to_base(250, 2, Optim_Info{}));

        for(size_t i = 0; i < random_runs; i++)
        {
            Optim_Info optims = generate_random_optims(generator);
            umax num = num_dist(*generator);
            Num base = base_dist(*generator);

            assert(auto_test_to_base(num, base, optims));
            assert(test_to_and_fro(num, base, optims));
        };
    }

    template <typename T>
    runtime_proc test_pow(Memory_Resource* memory, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        enum Pow_Algorhitm
        {
            OPTIMAL,
            TRIVIAL,
            SQUARING
        };

        runtime_proc pow_ = [&](CSlice num, size_t pow, Optim_Info optims, Pow_Algorhitm algorhitm) -> Vector {
            size_t required_size = required_pow_out_size(num, pow);
            size_t aux_size = required_pow_by_squaring_aux_size(num, pow)*100 + 100;
            Vector out = make_sized_vector<T>(required_size, resource);
            Vector aux = make_sized_vector<T>(aux_size, resource);

            Slice out_s = to_slice(&out);
            Slice aux_s = to_slice(&aux);

            Slice powed;
            if(algorhitm == TRIVIAL)
                powed = ::trivial_pow<T>(&out_s, &aux_s, num, pow, optims);
            else if(algorhitm == SQUARING)
                powed = ::pow_by_squaring<T>(&out_s, &aux_s, num, pow, optims);
            else
                powed = ::pow<T>(&out_s, &aux_s, num, pow, optims);

            resize(&out, powed.size);
            return out;
        };

        runtime_proc pow = [&](umax num, size_t pow, Optim_Info optims, Pow_Algorhitm algorhitm) -> umax {
            Vector num_ = make_vector_of_digits<T>(num, resource);
            Vector powed = pow_(to_slice(num_), pow, optims, algorhitm);
            CSlice powed_s = to_slice(powed);

            umax res = unwrap(to_number(powed_s));
            return res;
        };

        runtime_proc pow_by_squaring = [&](umax num, umax power, Optim_Info optims = Optim_Info{}) -> umax {
            return pow(num, power, optims, SQUARING);
        };

        runtime_proc trivial_pow = [&](umax num, umax power, Optim_Info optims = Optim_Info{}) -> umax {
            return pow(num, power, optims, TRIVIAL);
        };

        runtime_proc optimal_pow = [&](umax num, umax power, Optim_Info optims = Optim_Info{}) -> umax {
            return pow(num, power, optims, OPTIMAL);
        };

        runtime_proc test_pow = [&](CSlice num, size_t pow, CSlice expected, Pow_Algorhitm algorhitm, Optim_Info optims) -> bool {
            Vector powed = pow_(num, pow, optims, algorhitm);
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

        assert(optimal_pow(0, 0) == 1);
        assert(optimal_pow(0, 1) == 0);
        assert(optimal_pow(0, 54) == 0);
        assert(optimal_pow(0, 5554135) == 0);
        assert(optimal_pow(2, 0) == 1);
        assert(optimal_pow(1048576, 0) == 1);
        assert(optimal_pow(5465313, 1) == 5465313);
        assert(optimal_pow(5465313, 2) == 5465313ull*5465313ull);
        assert(optimal_pow(2, 10) == 1024);
        assert(optimal_pow(2, 20) == 1048576);
        assert(optimal_pow(3, 10) == 59049);
        assert(optimal_pow(3, 20) == 59049ull*59049ull);
        assert(optimal_pow(9, 10) == 59049ull*59049ull);
        assert(optimal_pow(5465, 3) == 5465ull*5465ull*5465ull);
    }
    

    template <typename T>
    runtime_proc test_root(Memory_Resource* memory, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using umax = umax;
        using Vector = Vector<T>;
        using Padded = Padded_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        Memory_Resource* resource = memory;

        runtime_proc root_ = [&](CSlice num, size_t root, Optim_Info optims) -> Vector {
            size_t required_size = root_required_out_size(log2(num), root, BIT_SIZE<T>);
            size_t aux_size = root_required_aux_size(log2(num), root, BIT_SIZE<T>);
            Vector out = make_sized_vector<T>(required_size, resource);
            Vector aux = make_sized_vector<T>(aux_size, resource);

            Slice out_s = to_slice(&out);
            Slice aux_s = to_slice(&aux);

            Slice rooted = ::root<T>(&out_s, &aux_s, num, root, optims);

            resize(&out, rooted.size);
            return out;
        };

        runtime_proc root = [&](umax num, size_t root, Optim_Info optims = Optim_Info{}) -> umax {
            Vector num_ = make_vector_of_digits<T>(num, resource);
            Vector rooted = root_(to_slice(num_), root, optims);
            CSlice rooted_s = to_slice(rooted);

            umax res = unwrap(to_number(rooted_s));
            return res;
        };


        assert(root(0, 0) == 1);
        assert(root(1, 0) == 1);
        assert(root(29, 3) == 3);
        assert(root(34, 5) == 2);
        assert(root(15625, 3) == 25);
        assert(root(15637, 3) == 25);
        assert(root(4096, 6) == 4);

        if(sizeof(T) != 1) //17th root is bigger than half bits of u8!
        {
            assert(root(18446744073709551614, 16) == 15);
            assert(root(18446744073709551614, 17) == 13);
        }
    }

    runtime_proc run_typed_tests(Memory_Resource* memory, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        std::cout << "TESTING UNTYPED "<< std::endl;
        test_misc(generator, random_runs);
        std::cout << "=== OK ===\n" << std::endl;
    }

    template <typename T>
    runtime_proc run_untyped_tests(Memory_Resource* memory, Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        std::cout << "TESTING: " << meta::type_name<T>() << std::endl;

        std::cout << "add...        ";
        test_add_overflow<T>(memory, generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "sub...        ";
        test_sub_overflow<T>(memory, generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "complement... ";
        test_complement_overflow<T>(memory, generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "shift...      ";
        test_shift_overflow<T>(memory, generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "mul short...  ";
        test_mul_overflow<T>(memory, generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;

        std::cout << "div short...  ";
        test_div_overflow_low<T>(memory, generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;

        std::cout << "mul long...   ";
        test_mul_quadratic<T>(memory, generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "mul long...   ";
        test_div_bit_by_bit<T>(memory, generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "pow...        ";
        test_pow<T>(memory, generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;

        std::cout << "root...       ";
        test_root<T>(memory, generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;

        std::cout << "conversion... ";
        const size_t quarter_runs = random_runs / 4;
        if constexpr(sizeof(T) >= sizeof(u8))
            test_to_base<T, u8>(memory, generator, quarter_runs);

        if constexpr(sizeof(T) >= sizeof(u16))
            test_to_base<T, u16>(memory, generator, quarter_runs);

        if constexpr(sizeof(T) >= sizeof(u32))
            test_to_base<T, u32>(memory, generator, quarter_runs);

        if constexpr(sizeof(T) >= sizeof(u64))
            test_to_base<T, u64>(memory, generator, quarter_runs);

        std::cout << "ok" << std::endl;
        std::cout << "=== OK ===\n" << std::endl;
    }


    void print(const char* name, u64 const& iters)
    {
        std::cout << name << std::endl;
        std::cout << "iters: " << iters << std::endl;
        std::cout << std::endl;
    }

    runtime_proc run_tests()
    {
        std::random_device os_seed;
        const u32 seed = os_seed();

        Random_Generator generator(seed);
        const size_t random_runs = 10000;
        const size_t controlled_runs = 500;

        auto ns = ellapsed_time([&](){
            #ifdef USE_CUSTOM_LIB
                jot::Allocator* def = jot::memory_globals::default_allocator();

                auto memory_res = def->allocate(jot::memory_constants::MEBI_BYTE*8, 8);
                defer(def->deallocate(memory_res.items, 8));

                jot::Ring_Allocator ring(memory_res.items, def);
                Memory_Resource* res = &ring;
            #else
                Memory_Resource* res = std::pmr::new_delete_resource();
            #endif

            run_typed_tests(res, &generator, random_runs, controlled_runs);
            run_untyped_tests<u8>(res, &generator, random_runs, controlled_runs);
            run_untyped_tests<u16>(res, &generator, random_runs, controlled_runs);
            run_untyped_tests<u32>(res, &generator, random_runs, controlled_runs);
            run_untyped_tests<u64>(res, &generator, random_runs, controlled_runs);
        });
        std::cout << "tests ellapsed: " << ns << " ns" << std::endl;
        
    }
}




#undef let 
#undef mut 
#undef proc 
#undef func 
#undef runtime_proc 
#undef runtime_func 
#undef cast