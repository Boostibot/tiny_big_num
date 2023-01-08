#pragma once

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <cctype>
#include <cmath>
#include <type_traits>

#include <iostream>
#include <random>
#include <limits>
#include <initializer_list>
#include <memory_resource>
#include <array>
#include <algorithm>

#include "pod_vector.h"
#include "benchmark.h"

#define TINY_NUM_TYPE_1 uint8_t
#define TINY_NUM_LIB_NAME tiny_num_u8
#include "tiny_big_num.h"

#define TINY_NUM_TYPE_2 uint16_t
#undef TINY_NUM_LIB_NAME
#define TINY_NUM_LIB_NAME tiny_num_u16
#include "tiny_big_num.h"

#define TINY_NUM_TYPE_3 uint32_t
#undef TINY_NUM_LIB_NAME
#define TINY_NUM_LIB_NAME tiny_num_u32
#include "tiny_big_num.h"

#define TINY_NUM_TYPE_4 uint64_t
#undef TINY_NUM_LIB_NAME
#define TINY_NUM_LIB_NAME tiny_num_u64
#include "tiny_big_num.h"

#define let const auto
#define mut auto
#define proc constexpr auto
#define func [[nodiscard]] constexpr auto
#define runtime_proc [[nodiscard]] auto
#define runtime_func auto

#define cast(...) (__VA_ARGS__)

namespace tiny_num_test
{
    using namespace helpers;

    namespace tiny_num
    {
        using namespace tiny_num_u8;
        using namespace tiny_num_u16;
        using namespace tiny_num_u32;
        using namespace tiny_num_u64;
    };

    using u8 = uint8_t;
    using u16 = uint16_t;
    using u32 = uint32_t;
    using u64 = uint64_t;

    using i8 = int8_t;
    using i16 = int16_t;
    using i32 = int32_t;
    using i64 = int64_t;

    using b8 = bool;
    using b16 = uint16_t;
    using b32 = uint32_t;
    using b64 = uint64_t;

    using f32 = float;
    using f64 = double;

    using cstring = const char*;
    using isize = i64;
    using usize = size_t;

    using umax = u64;

    template <typename T>
    constexpr size_t BIT_SIZE = sizeof(T) * CHAR_BIT;

    using std::move;
    using std::size;
    using std::data;

    template<typename T>
    func make_sized_vector(size_t size) -> POD_Vector<T>
    {
        POD_Vector<T> vec;
        vec.resize(size);
        return vec;
    }

    template<typename T>
    func push_back(POD_Vector<T>* vec, T pushed)
    {
        vec->push_back(move(pushed));
    }

    template<typename T>
    func resize(POD_Vector<T>* vec, size_t size)
    {
        vec->resize(size);
    }

    template<typename T>
    func slice(POD_Vector<T> const& vec) -> Slice<const T>
    {
        return Slice<const T>{data(vec), cast(usize) size(vec)};
    }

    func slice(cstring str) -> Slice<const char>
    {
        size_t size = 0;
        while(str[size] != '\0')
            size ++;
        
        return Slice<const char>{str, size};
    }

    template<typename T>
    func slice(POD_Vector<T>* vec) -> Slice<T>
    {
        return Slice<T>{std::data(*vec), cast(usize) size(*vec)};
    }

    template<typename T>
    func make_vector_of_slice(Slice<const T> const& items) -> POD_Vector<T>
    {
        mut vec = make_sized_vector<T>(items.size);
        mut vec_slice = slice(&vec);
        tiny_num::copy_slice(&vec_slice, items);

        return vec;
    }

    template<typename T>
    func make_vector_of_digits(umax val) -> POD_Vector<T>
    {
        constexpr size_t count = MAX_TYPE_SIZE_FRACTION<T>;
        POD_Vector<T> vec = make_sized_vector<T>(count);
        Slice<T> digits = slice(&vec);
        Slice<T> converted_slice = tiny_num::from_number(&digits, val);
        resize(&vec, converted_slice.size);
        return vec;
    }

    template <typename T>
    struct Padded_POD_Vector
    {
        POD_Vector<T> vector;
        size_t prefix_size = 0;
        size_t content_size = 0;
        size_t postfix_size = 0;
    };

    template <typename T>
    func make_padded_vector_of_digits(umax val, size_t pref_size = 1, size_t post_size = 1) -> Padded_POD_Vector<T>
    {
        constexpr size_t count = MAX_TYPE_SIZE_FRACTION<T>;
        size_t total_size = pref_size + count + post_size;

        Padded_POD_Vector<T> padded;

        padded.vector = make_sized_vector<T>(total_size);
        Slice<T> whole = slice(&padded.vector);
        Slice<T> digits = slice_size(whole, pref_size, count);
        Slice<T> content = tiny_num::from_number(&digits, val);

        padded.prefix_size = pref_size;
        padded.content_size = content.size;
        padded.postfix_size = post_size;

        return padded;
    }

    template <typename T>
    func content(Padded_POD_Vector<T>* padded) -> Slice<T>
    {
        Slice<T> whole = slice(&padded->vector);
        return slice_size(whole, padded->prefix_size, padded->content_size);
    }

    template <typename T>
    func content(Padded_POD_Vector<T> const& padded) -> Slice<const T>
    {
        Slice<const T> whole = slice(padded.vector);
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
            double val = std::pow(this->quotient, this->uniform_distribution(gen)) - this->shift;
            if constexpr (std::is_floating_point_v<T>)
                return cast(T) val;
            else
                return cast(T) round(val);
        }

    };

    template <typename T, typename ReturnT = umax>
    runtime_func make_exponential_distribution(ReturnT min = 0, ReturnT max = std::numeric_limits<ReturnT>::max()) -> Uniform_Exponential_Distribution<ReturnT>
    {
        double quot = std::pow(2, BIT_SIZE<T>);
        return {cast(double) min, cast(double) max, quot};
    }

    runtime_func generate_random_optims(Random_Generator* gen) -> Optims
    {
        let exp_dist = make_exponential_distribution<u8, size_t>(0, 100);
        //std::uniform_int_distribution<> dist{0, 100};
        std::uniform_int_distribution<> uniform_bool{0, 1};
        let random_bool = [&]() -> bool { return cast(bool) uniform_bool(*gen); };
        let random_num = [&]() -> size_t { return cast(size_t) exp_dist(*gen); };

        Optims optims;
        optims.mul_shift = random_bool();
        optims.mul_consts = random_bool();
        optims.mul_half_bits = random_bool();

        optims.div_shift = random_bool();
        optims.div_consts = random_bool();
        optims.rem_optims = random_bool();

        optims.max_recursion_depth = random_num();
        optims.mul_quadratic_both_below_size = random_num();
        optims.mul_quadratic_single_below_size = random_num();
        optims.pow_trivial_below_power = random_num();

        return optims;
    }

    runtime_proc test_misc(Random_Generator* generator, size_t random_runs)
    {
        {
            using namespace tiny_num_u64;

            assert(find_last_set_bit(0) == 0);
            assert(find_last_set_bit(1) == 0);
            assert(find_last_set_bit(2) == 1);
            assert(find_last_set_bit(3) == 1);
            assert(find_last_set_bit(4) == 2);
            assert(find_last_set_bit(7) == 2);
            assert(find_last_set_bit(8) == 3);
            assert(find_last_set_bit(76) == 6);
            assert(find_last_set_bit(513) == 9);
            assert(find_last_set_bit(641) == 9);
            assert(find_last_set_bit(1024) == 10);
            assert(find_last_set_bit(0xf4f1f63db9b7e800) == 63);
            assert(find_last_set_bit(0b010100) == 4);
            assert(find_last_set_bit(0b010000) == 4);

            assert(find_first_set_bit(0) == 0);
            assert(find_first_set_bit(0b000001) == 0);
            assert(find_first_set_bit(0b010100) == 2);
            assert(find_first_set_bit(0b010101) == 0);
            assert(find_first_set_bit(0b010000) == 4);
            assert(find_first_set_bit(0b010100) == 2);
            assert(find_first_set_bit(0b010010) == 1);
        
            assert(pop_count(0b0000'0000) == 0);
            assert(pop_count(0b0000'0001) == 1);
            assert(pop_count(0b1000'0000) == 1);
            assert(pop_count(0b1010'0000) == 2);
            assert(pop_count(0b1010'0011) == 4);
            assert(pop_count(0b0111'1111) == 7);
            assert(pop_count(cast(u32) -1) == 32);
            assert(pop_count(cast(umax) -1) == 64);

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
            assert(single_root_shifting(18446744073709551614, 17ull) == 13);

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

            assert((highest_power_of_in_bits(2, 2) == Power_And_Powed{1, 2}));
            assert((highest_power_of_in_bits(1, 2) == Power_And_Powed{0, 1}));
            assert((highest_power_of_in_bits(16, 2) == Power_And_Powed{15, 1 << 15}));
            assert((highest_power_of_in_bits(8, 3) == Power_And_Powed{5, 243}));

            assert((highest_power_of_in_bits(0, 10) == Power_And_Powed{0, 1}));
            assert((highest_power_of_in_bits(1, 10) == Power_And_Powed{0, 1}));
            assert((highest_power_of_in_bits(2, 10) == Power_And_Powed{0, 1}));
            assert((highest_power_of_in_bits(4, 10) == Power_And_Powed{1, 10}));
            assert((highest_power_of_in_bits(8, 10) == Power_And_Powed{2, 100}));
            assert((highest_power_of_in_bits(16, 10) == Power_And_Powed{4, 10'000}));
            assert((highest_power_of_in_bits(32, 10) == Power_And_Powed{9, 1'000'000'000}));
            assert((highest_power_of_in_bits(33, 10) == Power_And_Powed{9, 1'000'000'000}));
            assert((highest_power_of_in_bits(34, 10) == Power_And_Powed{10, 10'000'000'000}));
            assert((highest_power_of_in_bits(64, 10) == Power_And_Powed{19, 10'000'000'000'000'000'000}));

        }

        assert(tiny_num_u8::high_mask() == 0xF0);
        assert(tiny_num_u8::high_mask(2) == 0b1111'1100);
        assert(tiny_num_u8::high_mask(5) == 0b1110'0000);
        assert(tiny_num_u16::high_mask() == 0xFF00);
        assert(tiny_num_u32::high_mask() == 0xFFFF0000);
        assert(tiny_num_u64::high_mask() == 0xFFFFFFFF00000000);

        assert(tiny_num_u8::low_mask() == 0x0F);
        assert(tiny_num_u8::low_mask(2) == 0b0000'0011);
        assert(tiny_num_u8::low_mask(5) == 0b0001'1111);
        assert(tiny_num_u16::low_mask() == 0x00FF);
        assert(tiny_num_u32::low_mask() == 0x0000FFFF);
        assert(tiny_num_u64::low_mask() == 0x00000000FFFFFFFF);

        assert(tiny_num_u32::low_bits(0b1010'1100'0001, 1) == 1);
        assert(tiny_num_u32::low_bits(0b1010'1100'0001, 2) == 1);
        assert(tiny_num_u32::low_bits(0b1010'1100'0001, 4) == 1);
        assert(tiny_num_u32::low_bits(0b1010'1100'0001, 6) == 1);
        assert(tiny_num_u32::low_bits(0b1010'1100'0001, 8) == 0b1100'0001);
        assert(tiny_num_u32::low_bits(1, 16) == 1);
        assert(tiny_num_u32::low_bits(1, 7) == 1);
        assert(tiny_num_u32::low_bits(1, 8) == 1);

        constexpr u64 max = std::numeric_limits<u64>::max()/2;
        constexpr u64 min = 2;
        let val_dist = make_exponential_distribution<u64>(min, max);

        assert(cast(u64) 1 << 42 == 4398046511104);
        assert(tiny_num_u64::single_log_heyley(4398046511104ull, 2ull) == 42);

        for(size_t i = 0; i < random_runs; i++)
        {
            u64 val = cast(u64) val_dist(*generator);
            let power_dist = make_exponential_distribution<u64>(min, val);
            u64 power = cast(u64) power_dist(*generator);

            u64 rooted = tiny_num_u64::single_root_shifting(val, power);
            u64 powed_low = tiny_num_u64::single_pow_by_squaring(rooted, power);
            u64 back_down = tiny_num_u64::single_root_shifting(powed_low, power);
            if(rooted > 1 && powed_low < 43980465111) //so that no overflow
            {
                u64 loged_heyley = tiny_num_u64::single_log_heyley(powed_low, rooted);
                u64 loged_bsearch = tiny_num_u64::single_log_heyley(powed_low, rooted);
                assert(loged_heyley == power);
                assert(loged_bsearch == power);
            }

            if(power < 100)
            {
                u64 ref_powed_low = tiny_num_u64::single_pow_trivial(rooted, power);
                assert(ref_powed_low == powed_low);
            }
            else
            {
                //most of the time (maybe always) the diff is exactly zero 
                // yet we still want to be sure this test wont break randomly
                double epsilon = 0.00000000001;
                double ref_ipow = std::pow(rooted, power);
                double powed_low_d = cast(double) powed_low;

                double diff = abs(ref_ipow - powed_low_d);
                double frac = diff / cast(double) val; 
                //we divide by val because the only way errors can happen is double conversion on
                // large numbers => we need to only consider large numbers and their 'relative' error
                assert(frac < epsilon);
            }
            assert(powed_low <= val);
            assert(back_down == rooted);


            u64 rooted_again = tiny_num_u64::single_root_shifting(rooted, power);
            u64 rooted_again_newton = tiny_num_u64::single_root_newton(rooted, power);

            assert(rooted_again == rooted_again_newton);
        }
    }

    template <typename T>
    runtime_proc test_add_overflow(Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        runtime_proc test_add_overflow_batch = [&](umax left_, umax right_, umax carry_in, Res expected, bool check_overflow = true) -> bool {
            POD_Vector left = make_vector_of_digits<T>(left_);
            POD_Vector right = make_vector_of_digits<T>(right_);
            POD_Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

            Padded pad_left = make_padded_vector_of_digits<T>(left_, 1);
            Padded pad_right = make_padded_vector_of_digits<T>(right_, 1);

            Slice out_s = slice<T>(&out);
            Slice pad_left_s = slice<T>(&pad_left.vector);
            Slice pad_right_s = slice<T>(&pad_right.vector);

            let loc = Location::OUT_OF_PLACE;

            let res1 = batch_add_overflow_long(&out_s,  slice(left), slice(right), loc, cast(T) carry_in);
            let res2 = batch_add_overflow_long(&pad_left_s, content(pad_left), slice(right), loc, cast(T) carry_in); //test self assign
            let res3 = batch_add_overflow_long(&pad_right_s, slice(left), content(pad_right), loc, cast(T) carry_in); 

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

            let num1 = to_number(res1_s);
            let num2 = to_number(res2_s);
            let num3 = to_number(res3_s);

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

            mut left1 = make_vector_of_digits<T>(adjusted_left); 
            mut left2 = make_vector_of_digits<T>(adjusted_left);

            Slice left1_s = slice(&left1);
            Slice left2_s = slice(&left2);

            let res2 = batch_add_overflow_short(&left1_s, left1_s, adjusted_single, Location::IN_PLACE);
            let res1 = batch_add_overflow_long(&left2_s, left2_s, Slice{&adjusted_single, 1}, Location::IN_PLACE);

            bool slices = is_equal(res1.slice, res2.slice);
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
    runtime_proc test_sub_overflow(Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        runtime_proc test_sub_overflow_batch = [&](umax left_, umax right_, umax carry_in, Res expected, bool check_overflow = true) -> bool {
            POD_Vector left = make_vector_of_digits<T>(left_);
            POD_Vector right = make_vector_of_digits<T>(right_);
            POD_Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

            Padded pad_left = make_padded_vector_of_digits<T>(left_, 1);
            Padded pad_right = make_padded_vector_of_digits<T>(right_, 1);

            Slice out_s = slice<T>(&out);
            Slice pad_left_s = slice<T>(&pad_left.vector);

            let loc = Location::OUT_OF_PLACE;

            let res1 = batch_sub_overflow_long(&out_s,  slice(left), slice(right), loc, cast(T) carry_in);
            let res2 = batch_sub_overflow_long(&pad_left_s, content(pad_left), slice(right), loc, cast(T) carry_in); //test self assign

            Slice res1_s = res1.slice;
            Slice res2_s = res2.slice;

            if(check_overflow == false)
                assert(left_ >= right_ && "must not overflow when auto testing");

            let num1 = to_number(res1_s);
            let num2 = to_number(res2_s);

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
            umax adjusted_carry = carry_in > tiny_num_u64::low_mask();

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
    runtime_proc test_complement_overflow(Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        //we dont test oveflow of this op 
        runtime_proc batch_complement_overflow = [&](umax left_, umax carry_in = 1) -> umax{
            Padded pad_left = make_padded_vector_of_digits<T>(left_, 1);
            POD_Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

            Slice left_s = content(&pad_left);
            Slice out_s = slice<T>(&out);
            Slice pad_left_s = slice<T>(&pad_left.vector);

            let loc = Location::OUT_OF_PLACE;

            let res1 = tiny_num::batch_complement_overflow(&out_s, left_s, cast(T) carry_in);
            let res2 = tiny_num::batch_complement_overflow(&pad_left_s, left_s, cast(T) carry_in); //test self assign

            let num1 = to_number(res1.slice);
            let num2 = to_number(res2.slice);

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
    runtime_proc test_shift_overflow(Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        runtime_proc shift_both = [&](umax left_, umax right, umax carry_in, bool is_shift_up, Iter_Direction direction) -> Batch_Op_Result{
            Padded pad_left = make_padded_vector_of_digits<T>(left_, 1);
            POD_Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

            Slice left_s = content(&pad_left);
            Slice out_s = slice<T>(&out);
            Slice whole = slice<T>(&pad_left.vector);

            Slice pad_left_s;
            if(direction == Iter_Direction::FORWARD)
                pad_left_s = whole;
            else
                pad_left_s = slice(whole, 1);

            let res1 = (is_shift_up) 
                ? batch_shift_up_overflow(&out_s, left_s, right, direction, cast(T) carry_in)
                : batch_shift_down_overflow(&out_s, left_s, right, direction, cast(T) carry_in);
            let res2 = (is_shift_up) 
                ? batch_shift_up_overflow(&pad_left_s, left_s, right, direction, cast(T) carry_in)
                : batch_shift_down_overflow(&pad_left_s, left_s, right, direction, cast(T) carry_in);

            let num1 = to_number(res1.slice);
            let num2 = to_number(res2.slice);

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
        using namespace tiny_num_u64;
        // (by doing find_last_set_bit) and then shifting it so that they will not overflow

        umax adjusted_left = left;
        umax adjusted_right = right;

        size_t log_left = log2(left) + 1; 
        size_t log_right = log2(right) + 1; 
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

        assert(log2(adjusted_left) + log2(adjusted_right) + 2 <= BIT_SIZE<umax>);
        return {adjusted_left, adjusted_right};
    }

    template <typename T>
    runtime_proc test_mul_overflow(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        runtime_proc test_mul_overflow_batch = [&](umax left_, umax right, umax carry_in, Res expected, Optims optims, bool check_overflow = true) -> bool {
            Padded pad_left = make_padded_vector_of_digits<T>(left_, 1);
            POD_Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

            Slice left_s = content(&pad_left);
            Slice out_s = slice<T>(&out);
            Slice pad_left_s = slice<T>(&pad_left.vector);

            let res1 = batch_mul_overflow(&out_s, left_s, cast(T) right, optims, cast(T) carry_in);
            let res2 = batch_mul_overflow(&pad_left_s, left_s, cast(T) right, optims, cast(T) carry_in);

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

            let num1 = to_number(res1_s);
            let num2 = to_number(res2_s);

            bool outputs_match = 
                num1 == expected.value && 
                num2 == expected.value; 

            bool overflows_match = res1.overflow == res2.overflow; 
            bool sizes_match = res1.slice.size == res2.slice.size; 

            if(check_overflow)
                overflows_match = overflows_match && res1.overflow == expected.overflow;

            return outputs_match && overflows_match && sizes_match;
        };

        runtime_proc auto_test_mul = [&](umax left, umax right, Optims optims) -> bool
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

                assert(test_mul_overflow_batch(0xFFEE, 0x80, 0, Res{0xF700, 0x7F}, optims));
                assert(test_mul_overflow_batch(0xFFEE, 0x80, 0, Res{0xF700, 0x7F}, optims));
                assert(test_mul_overflow_batch(0xA1, 0x11, 0, Res{0xB1, 0xA}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 1, 0, Res{0xABCDEF00, 0}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 0, 0, Res{0, 0}, optims));
                assert(test_mul_overflow_batch(0xABCDEF00, 0xAB, 0, Res{0xC28EA500, 0x72}, optims));
                assert(test_mul_overflow_batch(13745, 137, 0, Res{0xBBB9, 0x1C}, optims));
            }
        }

        assert(auto_test_mul(0x14156f, 0x3c, Optims{}));

        let exponential = make_exponential_distribution<T, umax>();
        let uniform = std::uniform_int_distribution<long long>(0, T(-1) >> 1);
        for(size_t i = 0; i < random_runs; i++)
        {
            umax left = exponential(*generator);
            umax right = uniform(*generator);
            let optims = generate_random_optims(generator);
            assert(auto_test_mul(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_div_overflow(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        //assert((single_div_overflow(cast(T) 0x12, 0x0A, 0) == Single_Overflow{1, 0x08}));

        using namespace tiny_num;

        runtime_proc test_div_overflow_batch = [&](umax left_, umax right, umax carry_in, Res expected, Optims optims) -> bool {
            constexpr size_t padding_before = 2;

            Padded pad_left = make_padded_vector_of_digits<T>(left_, 1);
            POD_Vector out = make_sized_vector<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

            Slice left_s = content(&pad_left);
            Slice out_s = slice<T>(&out);
            Slice pad_left_s = slice(slice<T>(&pad_left.vector), 1);

            let res1 = batch_div_overflow(&out_s, left_s, cast(T) right, optims, cast(T) carry_in);
            let res2 = batch_div_overflow(&pad_left_s, left_s, cast(T) right, optims, cast(T) carry_in);

            let num1 = to_number(res1.slice);
            let num2 = to_number(res2.slice);

            bool outputs_match = 
                num1 == expected.value && 
                num2 == expected.value; 

            bool overflows_match = 
                res1.overflow == expected.overflow && 
                res2.overflow == expected.overflow; 

            bool sizes_match = res1.slice.size == res2.slice.size; 

            return outputs_match && overflows_match && sizes_match;
        };
    
        runtime_proc auto_test_div = [&](umax left, umax right, Optims optims) -> bool {
            if(right == 0)
                return true;
        
            umax adjusted_right = tiny_num_u64::low_bits(cast(T) right);
            umax res = left / adjusted_right;
            umax rem = left % adjusted_right;

            //carry in is basically added as an extra digit to the divided number
            // its entire purpose is linking multiple div blocks if desired (ie we store number chunked for some reason)
            // we dont test this yet because the setup would basically exactly mirror of the operation
            // (because its dependent on the template type)

            return test_div_overflow_batch(left, adjusted_right, 0, Res{res, rem}, optims);
        };

        if constexpr(std::is_same_v<T, u8>)
        {
            for(size_t i = 0; i < controlled_runs; i++)
            {
                let optims = generate_random_optims(generator);
                assert(test_div_overflow_batch(0, 1, 0, Res{0, 0}, optims));
                assert(test_div_overflow_batch(1, 1, 0, Res{1, 0}, optims));
                assert(test_div_overflow_batch(25, 1, 0, Res{25, 0}, optims));
                assert(test_div_overflow_batch(25, 5, 0, Res{5, 0}, optims));
                assert(test_div_overflow_batch(0xFF, 0x08, 0, Res{31, 7}, optims));
                assert(test_div_overflow_batch(0xA1, 0x0F, 0, Res{10, 11}, optims));
                assert(test_div_overflow_batch(0xA1C0, 0x08, 0x1, Res{0x3438, 0}, optims));
                assert(test_div_overflow_batch(0xABCDEF00, 1, 0, Res{0xABCDEF00, 0}, optims));
                assert(test_div_overflow_batch(0x99AA, 0x08, 0x9, Res{0x3335, 2}, optims));
                assert(test_div_overflow_batch(0xABCDEF00, 0x07, 0xC, Res{0xcf668fdb, 0x3}, optims));
                //for carry in the result is obtained as follows:
                // left:  0x0ABCDEF00 
                // carry: 0xC00000000
                // eff:   0xCABCDEF00 
                // eff / right = 0x1CF668FDB => 0xCF668FDB (removal of the added - carried digit)
                // eff % right = 0x3
                assert(test_div_overflow_batch(0xABCDEF00, 0x07, 0, Res{411771428, 4}, optims));
                assert(test_div_overflow_batch(0xABCDEF00, 0x0F, 0, Res{192160000, 0}, optims));
                assert(test_div_overflow_batch(13745, 0xB, 0, Res{1249, 6}, optims));
            }
        }

        let exponential = make_exponential_distribution<T, umax>();
        let uniform = std::uniform_int_distribution<long long>(0, T(-1) >> (BIT_SIZE<T> / 2));
        for(size_t i = 0; i < random_runs; i++)
        {
            umax left = exponential(*generator);
            umax right = uniform(*generator);
            let optims = generate_random_optims(generator);
            assert(auto_test_div(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_mul(Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        runtime_proc test_fused_mul_add = [&](umax left_, umax right_, umax coef, umax expected, Optims optims = Optims{}) -> bool{
            POD_Vector left = make_vector_of_digits<T>(left_);
            POD_Vector right = make_vector_of_digits<T>(right_);
            POD_Vector output = make_sized_vector<T>(size(left));

            Slice left_s = slice(&left);
            Slice right_s = slice(&right);
            Slice output_s = slice(&output);

            let res_out = batch_fused_mul_add_overflow(&output_s, left_s, right_s, cast(T) coef, optims, Location::OUT_OF_PLACE);
            let res_in = batch_fused_mul_add_overflow(&left_s, left_s, right_s, cast(T) coef, optims, Location::IN_PLACE);
        
            push_back(&output, res_out.overflow);
            push_back(&left, res_in.overflow);

            left_s = slice(&left);
            output_s = slice(&output);

            let num_out = to_number(left_s);
            let num_in = to_number(output_s);

            return expected == num_out && expected == num_in; 
        };

        runtime_proc auto_test_fused_mul_add = [&](umax left, umax right, umax coef, Optims optims = Optims{}) -> bool{
            return test_fused_mul_add(left, right, coef, left + coef*right, optims);
        };

        runtime_proc test_mul = [&](umax left_, umax right_, umax expected, Optims optims = Optims{}) -> bool{
            POD_Vector left = make_vector_of_digits<T>(left_);
            POD_Vector right = make_vector_of_digits<T>(right_);

            size_t to_size = tiny_num_u64::required_mul_out_size(size(left), size(right));
            size_t aux_size1 = tiny_num_u64::required_mul_quadratic_aux_size(size(left), size(right));
            size_t aux_size2 = tiny_num_u64::required_mul_karatsuba_aux_size(size(left), size(right), 1);

            size_t aux_size = max(aux_size1, aux_size2)*10; //to allow all karatsuba recursions

            POD_Vector res = make_sized_vector<T>(to_size);
            POD_Vector aux = make_sized_vector<T>(aux_size);

            Slice res_s = slice(&res);
            Slice aux_s = slice(&aux);
            CSlice left_s = slice(left);
            CSlice right_s = slice(right);

            tiny_num::mul_quadratic(&res_s, &aux_s, left_s, right_s, optims);
            let normal = to_number(res_s);

            tiny_num::mul_quadratic_fused(&res_s, left_s, right_s, optims);
            let fused = to_number(res_s);

            tiny_num::mul_karatsuba(&res_s, &aux_s, left_s, right_s, optims);
            let karatsuba = to_number(res_s);

            return expected == normal && expected == fused && expected == karatsuba;
        };

        runtime_proc auto_test_mul = [&](umax left, umax right, Optims optims) -> bool
        {
            let adjusted = adjust_to_not_mul_overflow(left, right);
            return test_mul(adjusted.left, adjusted.right, adjusted.left * adjusted.right, optims);
        };

        Optims shift = {true, true, false};

        assert(test_fused_mul_add(1, 1, 1, 2));
        assert(test_fused_mul_add(1, 1, 1, 2, shift));
        assert(test_fused_mul_add(1, 1, 5, 6));
        assert(test_fused_mul_add(1, 1, 5, 6, shift));
        assert(test_fused_mul_add(10, 1, 5, 15, shift));
        assert(auto_test_fused_mul_add(0xFFFF, 1, 0x80, shift));
        assert(auto_test_fused_mul_add(0x1000000000, 0x574ce80, 0x10, shift));
        assert(auto_test_fused_mul_add(0x0da8824af2592c00, 0x0da8824af2592c00, 0x2, shift));
        
        assert(test_mul(0, 0, 0));
        assert(test_mul(1, 0, 0));
        assert(test_mul(0, 1, 0));
        assert(test_mul(1, 1, 1));
        assert(test_mul(25, 1, 25));
        assert(test_mul(25, 5, 125));
        assert(test_mul(0xFF, 0x10, 0xFF0));
        assert(test_mul(0xA1, 0x11, 0xAB1));
        assert(test_mul(0xABCDEF00, 1, 0xABCDEF00));
        assert(test_mul(0xABCDEF00, 0, 0));
        assert(test_mul(13745, 137, 0x1CBBB9));
        assert(test_mul(0xFFFF, 0x10000, 0xFFFF0000));
        assert(test_mul(0xABCD, 0x10001, 0xABCDABCD));
        assert(test_mul(0xAB, 0x1254, 0xAB*0x1254));
        assert(test_mul(0x1000, 0x1200, 0x1000*0x1200));
        assert(test_mul(0xABCD, 0x54, 0xABCD*0x54));
        assert(test_mul(0xABCD, 0x1254, 0xABCD*0x1254));
        assert(test_mul(0xABCDEF124, 0x1254, 0xC4CDA61BA7D0));
        assert(test_mul(0x1254, 0xABCDEF124, 0xC4CDA61BA7D0));

        assert(auto_test_mul(0x1d15a574ce80, 0x38, shift));
        assert(auto_test_mul(0x1d15a574ce80, 0xFFFF, shift));
        assert(auto_test_mul(0x1d15a574ce80, 0x1010, shift));
        assert(auto_test_mul(0x1d15a574ce80, 0x7838, shift));
        assert(auto_test_mul(0x1d15a574ce80, 0x27838, shift));
        assert(auto_test_mul(60136640465, 262154, shift));

        let distribution = make_exponential_distribution<T>();
        for(size_t i = 0; i < random_runs; i++)
        {
            Optims optims = generate_random_optims(generator);
            umax left = distribution(*generator);
            umax right = distribution(*generator);

            assert(auto_test_mul(left, right, optims));
        };
    }

    template <typename T>
    runtime_proc test_div(Random_Generator* generator, size_t random_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        runtime_proc test_div = [&](umax num_, umax den_, umax ex_quotient, umax ex_remainder, Optims optims = Optims{}) -> bool{
            POD_Vector num = make_vector_of_digits<T>(num_); 
            POD_Vector den = make_vector_of_digits<T>(den_);
            size_t quo_size = tiny_num_u64::required_div_quotient_size(size(num), size(den));
            size_t rem_size = tiny_num_u64::required_div_remainder_size(size(num), size(den));

            POD_Vector quo = make_sized_vector<T>(quo_size); 
            POD_Vector rem = make_sized_vector<T>(rem_size);

            Slice quo_s = slice(&quo);
            Slice rem_s = slice(&rem);
            Slice num_s = slice(&num);
            Slice den_s = slice(&den);

            let div_res = tiny_num::div(&quo_s, &rem_s, num_s, den_s, optims);
            let quotient = to_number(div_res.quotient);
            let remainder = to_number(div_res.remainder);
            return ex_quotient == quotient && ex_remainder == remainder;
        };

        runtime_proc auto_test_div = [&](umax left, umax right, Optims optims = Optims{}) -> bool{
            if(right == 0)
                return true;;

            umax quotient = left / right;
            umax remainder = left % right;
            return test_div(left, right, quotient, remainder, optims);
        };

        assert(test_div(0, 1, 0, 0));
        assert(test_div(1, 1, 1, 0));
        assert(test_div(25, 1, 25, 0));
        assert(test_div(1, 25, 0, 1));
        assert(test_div(25, 26, 0, 25));
        assert(test_div(25, 5, 5, 0));
        assert(test_div(28, 5, 5, 3));
        assert(test_div(0xFF, 0x10, 0xF, 0xF));
        assert(test_div(183, 27, 6, 21));
        assert(test_div(0xABCDEF00, 1, 0xABCDEF00, 0));
        assert(test_div(0xABCDEF00, 0xABCDEF00, 1, 0));
        assert(test_div(0xABCDEF00, 0xABCDEF10, 0, 0xABCDEF00));
        assert(test_div(0xAB, 0xABCD, 0, 0xAB));
        assert(test_div(0xABCE, 0xABCD, 1, 1));

        assert(auto_test_div(13745, 6435));
        assert(auto_test_div(4643477, 4796));
        assert(auto_test_div(7893, 01047));
        assert(auto_test_div(65536, 256));
        assert(auto_test_div(256, 65536));
        assert(auto_test_div(256*413541, 256));
        assert(auto_test_div(959464934, 1111));
        assert(auto_test_div(27417318639, 57432413));
        assert(auto_test_div(6994935563924682752, 67877739413342216));

        let distribution = make_exponential_distribution<T>();
        for(size_t i = 0; i < random_runs; i++)
        {
            Optims optims = generate_random_optims(generator);
            umax left = distribution(*generator);
            umax right = distribution(*generator);
            assert(auto_test_div(left, right, optims));
        };
    }

    template <typename T>
    func to_base(Slice<const T> num, T base, Optims optims) -> POD_Vector<char> {
        POD_Vector<T> temp = make_vector_of_slice<T>(num);
        size_t required_size = tiny_num_u64::required_to_base_out_size(num.size, base, BIT_SIZE<T>);
        POD_Vector<char> out = make_sized_vector<char>(required_size);

        Slice<T> temp_s = slice(&temp);
        Slice<char> out_s = slice(&out);

        Slice<char> converted = (tiny_num::to_chars(&out_s, &temp_s, num, base, optims));

        resize(&out, converted.size);
        return out;
    };

    template <typename T>
    func from_base(Slice<const char> rep, T base, Optims optims) -> POD_Vector<T> {
        size_t required_size = tiny_num_u64::required_from_base_out_size(rep.size, base, BIT_SIZE<T>);
        POD_Vector<T> out = make_sized_vector<T>(required_size + 1);

        Slice<T> out_s = slice(&out);
        let converted = tiny_num::from_chars(&out_s, rep, base, optims);
        assert(converted.ok && tiny_num::is_striped_number(converted.slice));
        resize(&out, converted.slice.size);
        return out;
    };

    template <typename T>
    func native_to_base(umax num, T base) -> POD_Vector<char> 
    {
        constexpr size_t num_digits = MAX_TYPE_SIZE_FRACTION<T>;
        const size_t max_size = tiny_num_u64::required_to_base_out_size(num_digits, base, BIT_SIZE<T>);

        POD_Vector<char> converted = make_sized_vector<char>(max_size);
        size_t size = 0;
        for(umax curr_val = num; curr_val != 0; curr_val /= base, size ++)
        {
            umax rem = curr_val % base;
            converted[size] = tiny_num_u64::to_char(rem);
        }

        resize(&converted, size);
        size_t half_size = size / 2;
        for(size_t i = 0; i < half_size; i++)
            std::swap(converted[i], converted[size - i - 1]);

        return converted;
    };

    template <typename T>
    func native_from_base(Slice<const char> rep, T base) -> umax 
    {
        umax value = 0;
        for(size_t i = 0; i < rep.size; i++)
        {
            umax muled = value * base;
            umax added = muled + tiny_num_u64::to_digit(rep[i]);
            assert(muled >= value);
            assert(added >= muled);

            value = added;
        }

        return value;
    };


    template <typename T>
    func slice_equals(Slice<T> left, Slice<T> right) -> bool
    {
        if(left.size != right.size)
            return false;

        for(size_t i = 0; i < left.size; i++)
            if(left[i] != right[i])
                return false;

        return true;
    }

    template <typename T>
    runtime_proc test_to_base(Random_Generator* generator, size_t random_runs)
    {

        runtime_proc manual_test_to_base = [&](umax num, T base, cstring expected_rep, Optims optims = Optims{}) -> bool 
        {
            POD_Vector<T> num_ = make_vector_of_digits<T>(num);
            POD_Vector<char> converted = to_base<T>(slice(num_), base, optims);

            return slice_equals(slice(converted), slice(expected_rep));
        };

        runtime_proc auto_test_to_base = [&](umax num, T base, Optims optims) -> bool 
        {
            POD_Vector<T> num_ = make_vector_of_digits<T>(num);
            POD_Vector<char> expected_rep = native_to_base<T>(num, base);
            POD_Vector<char> converted = to_base<T>(slice(num_), base, optims);

            return slice_equals(slice(converted), slice(expected_rep));
        };

        runtime_proc manual_test_from_base = [&](cstring rep, T base, umax expected_num, Optims optims = Optims{}) -> bool 
        {
            POD_Vector<T> expected_num_ = make_vector_of_digits<T>(expected_num);
            POD_Vector<T> converted = from_base<T>(slice(rep), base, optims);

            return slice_equals(slice(converted), slice(expected_num_));
        };

        runtime_proc auto_test_from_base = [&](Slice<const char> rep, T base, Optims optims) -> bool 
        {
            umax num = native_from_base<T>(rep, base);
            POD_Vector<T> expected_num = make_vector_of_digits<T>(num);
            POD_Vector<T> converted = from_base<T>(rep, base, optims);

            return slice_equals(slice(converted), slice(expected_num));
        };

        runtime_proc test_to_and_fro = [&](umax value, T base, Optims optims) -> bool
        {
            POD_Vector<T> initial = make_vector_of_digits<T>(value);
            let rep = to_base<T>(slice(initial), base, optims);
            let num = from_base<T>(slice(rep), base, optims);
            Slice<const T> num_s = slice(num);
            let obtained_value = tiny_num::to_number(num_s);

            bool match = tiny_num::is_equal(slice(num), slice(initial));
            return match;
        };

        assert((manual_test_to_base(1, 2, "1")));
        assert((manual_test_to_base(2, 2, "10")));
        assert((manual_test_to_base(4, 2, "100")));
        assert((manual_test_to_base(5, 2, "101")));
        assert((manual_test_to_base(10, 2, "1010")));

        assert((manual_test_to_base(0, 9, "")));
        assert((manual_test_to_base(1, 9, "1")));
        assert((manual_test_to_base(2, 9, "2")));
        assert((manual_test_to_base(65453, 9, "108705")));
        assert((manual_test_to_base(844354, 9, "1526211")));
        assert((manual_test_to_base(844354, 10, "844354")));

        assert((manual_test_from_base("1", 2, 1)));
        assert((manual_test_from_base("10", 2, 2)));
        assert((manual_test_from_base("100", 2, 4)));
        assert((manual_test_from_base("101", 2, 5)));
        assert((manual_test_from_base("1010", 2, 10)));

        assert((manual_test_from_base("", 9, 0)));
        assert((manual_test_from_base("1", 9, 1)));
        assert((manual_test_from_base("2", 9, 2)));
        assert((manual_test_from_base("108705", 9, 65453)));
        assert((manual_test_from_base("1526211", 9, 844354)));
        assert((manual_test_from_base("844354", 10, 844354)));

        assert((test_to_and_fro(10, 2, Optims{})));
        assert((test_to_and_fro(59180470887883776, 3, Optims{})));

        let num_dist = make_exponential_distribution<T>();

        //algorhitm requires the base to be low_bits only
        // this is only ever a problem for u8 (which is testing only)
        T max_base = sizeof(T) == 1 ? 15 : 36;

        let base_dist = make_exponential_distribution<T, T>(2, max_base);
        assert(auto_test_to_base(250, 2, Optims{}));

        for(size_t i = 0; i < random_runs; i++)
        {
            Optims optims = generate_random_optims(generator);
            umax num = num_dist(*generator);
            T base = base_dist(*generator);

            assert(auto_test_to_base(num, base, optims));
            assert(test_to_and_fro(num, base, optims));
        };
    }

    template <typename T>
    runtime_proc test_pow(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        enum Pow_Algorhitm
        {
            OPTIMAL,
            TRIVIAL,
            SQUARING
        };

        using namespace tiny_num;

        runtime_proc pow_ = [&](CSlice num, size_t pow, Optims optims, Pow_Algorhitm algorhitm) -> POD_Vector {
            size_t log2num = log2(num);

            size_t required_size = tiny_num_u64::required_pow_out_size(log2num, pow, BIT_SIZE<T>);
            size_t aux_size = tiny_num_u64::required_pow_by_squaring_aux_size(log2num, pow, BIT_SIZE<T>)*100 + 100;
            POD_Vector out = make_sized_vector<T>(required_size);
            POD_Vector aux = make_sized_vector<T>(aux_size);

            Slice out_s = slice(&out);
            Slice aux_s = slice(&aux);

            Slice powed;
            if(algorhitm == TRIVIAL)
                powed = tiny_num::pow_trivial(&out_s, &aux_s, num, pow, optims);
            else if(algorhitm == SQUARING)
                powed = tiny_num::pow_by_squaring(&out_s, &aux_s, num, pow, optims);
            else
                powed = tiny_num::pow(&out_s, &aux_s, num, pow, optims);

            resize(&out, powed.size);
            return out;
        };

        runtime_proc pow = [&](umax num, size_t pow, Optims optims, Pow_Algorhitm algorhitm) -> umax {
            POD_Vector num_ = make_vector_of_digits<T>(num);
            POD_Vector powed = pow_(slice(num_), pow, optims, algorhitm);
            CSlice powed_s = slice(powed);

            umax res = to_number(powed_s);
            return res;
        };

        runtime_proc pow_by_squaring = [&](umax num, umax power, Optims optims = Optims{}) -> umax {
            return pow(num, power, optims, SQUARING);
        };

        runtime_proc pow_trivial = [&](umax num, umax power, Optims optims = Optims{}) -> umax {
            return pow(num, power, optims, TRIVIAL);
        };

        runtime_proc optimal_pow = [&](umax num, umax power, Optims optims = Optims{}) -> umax {
            return pow(num, power, optims, OPTIMAL);
        };

        runtime_proc test_pow = [&](CSlice num, size_t pow, CSlice expected, Pow_Algorhitm algorhitm, Optims optims) -> bool {
            POD_Vector powed = pow_(num, pow, optims, algorhitm);
            CSlice powed_s = slice(powed);

            assert(is_striped_number(powed_s));
            assert(is_striped_number(expected));
        
            return is_equal(powed_s, expected);
        };

        assert(pow_trivial(0, 0) == 1);
        assert(pow_trivial(0, 1) == 0);
        assert(pow_trivial(0, 54) == 0);
        assert(pow_trivial(0, 5554135) == 0);
        assert(pow_trivial(2, 0) == 1);
        assert(pow_trivial(1048576, 0) == 1);
        assert(pow_trivial(5465313, 1) == 5465313);
        assert(pow_trivial(5465313, 2) == 5465313ull*5465313ull);
        assert(pow_trivial(5465, 3) == 5465ull*5465ull*5465ull);
        assert(pow_trivial(2, 10) == 1024);
        assert(pow_trivial(2, 20) == 1048576);

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
    runtime_proc test_root(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        using Res = Batch_Op_Result;
        using POD_Vector = POD_Vector<T>;
        using Padded = Padded_POD_Vector<T>;
        using CSlice = Slice<const T>;
        using Slice = Slice<T>;

        using namespace tiny_num;

        runtime_proc root_ = [&](CSlice num, size_t root, Optims optims) -> POD_Vector {
            size_t required_size = tiny_num_u64::required_root_out_size(log2(num), root, BIT_SIZE<T>);
            size_t aux_size = tiny_num_u64::required_root_aux_size(log2(num), root, BIT_SIZE<T>);
            POD_Vector out = make_sized_vector<T>(required_size);
            POD_Vector aux = make_sized_vector<T>(aux_size);

            Slice out_s = slice(&out);
            Slice aux_s = slice(&aux);

            Slice rooted = tiny_num::root(&out_s, &aux_s, num, root, optims);

            resize(&out, rooted.size);
            return out;
        };

        runtime_proc root = [&](umax num, size_t root, Optims optims = Optims{}) -> umax {
            POD_Vector num_ = make_vector_of_digits<T>(num);
            POD_Vector rooted = root_(slice(num_), root, optims);
            CSlice rooted_s = slice(rooted);

            umax res = to_number(rooted_s);
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
    runtime_proc run_typed_tests(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
    {
        std::cout << "TESTING UNTYPED "<< std::endl;
        test_misc(generator, random_runs);
        std::cout << "=== OK ===\n" << std::endl;
    }

    template <typename T>
    runtime_proc run_untyped_tests(Random_Generator* generator, size_t random_runs, size_t controlled_runs, cstring type_name)
    {
        std::cout << "TESTING: " << type_name << std::endl;


        std::cout << "add...        ";
        test_add_overflow<T>(generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "sub...        ";
        test_sub_overflow<T>(generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "complement... ";
        test_complement_overflow<T>(generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "shift...      ";
        test_shift_overflow<T>(generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "mul short...  ";
        test_mul_overflow<T>(generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;

        std::cout << "div short...  ";
        test_div_overflow<T>(generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;
        
        std::cout << "mul long...   ";
        test_mul<T>(generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "div long...   ";
        test_div<T>(generator, random_runs);
        std::cout << "ok" << std::endl;

        std::cout << "pow...        ";
        test_pow<T>(generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;

        std::cout << "root...       ";
        test_root<T>(generator, random_runs, controlled_runs);
        std::cout << "ok" << std::endl;

        std::cout << "conversion... ";
        test_to_base<T>(generator, random_runs);
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

        auto ns = benchmark::ellapsed_time([&](){
            run_typed_tests(&generator, random_runs, controlled_runs);
            run_untyped_tests<u8>(&generator, random_runs, controlled_runs, "u8");
            run_untyped_tests<u16>(&generator, random_runs, controlled_runs, "u16");
            run_untyped_tests<u32>(&generator, random_runs, controlled_runs, "u32");
            run_untyped_tests<u64>(&generator, random_runs, controlled_runs, "u64");
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
