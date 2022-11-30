#pragma once


#include <iostream>
#include <cctype>
#include <random>
#include <limits>
#include <initializer_list>
#include <algorithm>

#include "benchmark.h"
#include "big_int.h"
#include "big_int_ops.h"

using Random_Generator = std::mt19937;

struct Batch_Op_Result
{
    Max_Unsigned_Type value;
    Max_Unsigned_Type overflow;

    bool constexpr operator ==(Batch_Op_Result const&) const noexcept = default;
};

template<integral T>
func make_big_int_with_size(size_t out_size) -> Big_Int_<T>
{
    Big_Int_<T> out;
    resize(&out, out_size);
    return out;
}

template <typename T>
constexpr size_t MAX_TYPE_SIZE_FRACTION = sizeof(Max_Unsigned_Type) / sizeof(T);

template <typename T>
struct Uniform_Exponential_Distribution
{
    double shift;
    double quotient;
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


template <integral T, integral Of = Max_Unsigned_Type>
runtime_func make_max_distribution(Of min = 0, Of max = std::numeric_limits<Of>::max()) -> Uniform_Exponential_Distribution<Of>
{
    double quot = pow(2, BIT_SIZE<T>);
    return {cast(double) min, cast(double) max, quot};
}

runtime_func generate_random_optims(Random_Generator* gen) -> Optim_Info
{
    std::uniform_int_distribution<> dist{0, 1};
    let random_bool = [&](){ return cast(bool) dist(*gen); };

    Optim_Info info;
    info.mul_shift = random_bool();
    info.mul_consts = random_bool();
    info.mul_half_bits = random_bool();

    info.div_shift = random_bool();
    info.div_consts = random_bool();
    info.rem_optims = random_bool();
    return info;
}

template <typename T>
struct Padded_Big_Int
{
    Big_Int_<T> big_int;
    Slice<T> num;
    Slice<T> whole_num;
    size_t before_digits;
    size_t after_digits;
};

template <typename T>
func make_padded_big_int(Max_Unsigned_Type val, size_t to_size, size_t before_digits, size_t after_digits) -> Padded_Big_Int<T>
{
    Big_Int_<T> big;
    size_t raw_digits = max(digits_to_represent<T, Max_Unsigned_Type>(), to_size);
    resize(&big, raw_digits + before_digits + after_digits);

    Slice<T> whole_num = big;
    Slice<T> from_start = slice(whole_num, before_digits);
    Slice<T> converted = from_number(&from_start, val);

    return Padded_Big_Int<T>{
        std::move(big),
        converted,
        whole_num,
        before_digits, 
        after_digits
    };
}

template <typename T>
func make_padded_big_int(Max_Unsigned_Type val, size_t before_digits = 1, size_t after_digits = 1) -> Padded_Big_Int<T>
{
    return make_padded_big_int<T>(val, 0, before_digits, after_digits);
}

proc test_misc()
{
    using Max = Max_Unsigned_Type;

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
}

template <typename T>
runtime_proc test_add_overflow(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    proc test_add_overflow_batch = [](Max left_, Max right_, Max carry_in, Res expected, bool check_overflow = true) -> bool {
        Big left = left_;
        Big right = right_;

        mut pad_left = make_padded_big_int<T>(left_);
        mut pad_right = make_padded_big_int<T>(right_);
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

        Slice<T> out_s = out;
        Slice<T> pad_left_s = pad_left.whole_num;
        Slice<T> pad_right_s = pad_right.whole_num;
        let loc = Op_Location::OUT_OF_PLACE;

        let res1 = ::add_overflow_batch<T>(&out_s, left, right, loc, cast(T) carry_in);
        let res2 = ::add_overflow_batch<T>(&pad_left_s, pad_left.num, right, loc, cast(T) carry_in); //test self assign
        let res3 = ::add_overflow_batch<T>(&pad_right_s, left, pad_right.num, loc, cast(T) carry_in); 

        Slice<T> res1_s = res1.slice;
        Slice<T> res2_s = res2.slice;
        Slice<T> res3_s = res3.slice;

        if(check_overflow == false)
        {
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
    
    proc auto_test_add = [](Max left, Max right, Max carry_in = 0) -> bool {
        Max highest_one = cast(Max) 1 << (BIT_SIZE<Max> - 1);
        Max adjusted_left = left & ~highest_one;
        Max adjusted_right = right & ~highest_one;
        Max adjusted_carry = cast(Max) low_bits(cast(T) carry_in);
        Max expected = add_no_overflow(adjusted_left, adjusted_right, adjusted_carry);

        return test_add_overflow_batch(adjusted_left, adjusted_right, adjusted_carry, Res{expected, 0}, false);
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

    let distribution = make_max_distribution<T>();
    for(size_t i = 0; i < random_runs; i++)
    {
        Max left = distribution(*generator);
        Max right = distribution(*generator);
        Max carry = distribution(*generator);
        assert(auto_test_add(left, right, carry));
    };
}

template <typename T>
runtime_proc test_sub_overflow(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    proc test_sub_overflow_batch = [](Max left_, Max right_, Max carry_in, Res expected, bool check_overflow = true) -> bool {
        Big left = left_;
        Big right = right_;

        mut pad_left = make_padded_big_int<T>(left_);
        mut pad_right = make_padded_big_int<T>(right_);
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

        Slice<T> out_s = out;
        Slice<T> pad_left_s = pad_left.whole_num;
        let loc = Op_Location::OUT_OF_PLACE;

        let res1 = ::sub_overflow_batch<T>(&out_s, left, right, loc, cast(T) carry_in);
        let res2 = ::sub_overflow_batch<T>(&pad_left_s, pad_left.num, right, loc, cast(T) carry_in); //test self assign

        Slice<T> res1_s = res1.slice;
        Slice<T> res2_s = res2.slice;

        if(check_overflow == false)
            assert(left_ >= right_ && "must not overflow when auto testing");

        int x = 0x04b94100065dcc86 - 0x04b94100065dcc76;
        let num1 = unwrap(to_number(res1_s));
        let num2 = unwrap(to_number(res2_s));

        bool outputs_match =    num1 == expected.value && num2 == expected.value;
        bool overflows_match =  res1.overflow == res2.overflow;
        bool sizes_match =      res1.slice.size == res2.slice.size;

        if(check_overflow)
            overflows_match = overflows_match && res1.overflow == expected.overflow;

        return outputs_match && overflows_match && sizes_match;
    };

    proc auto_test_sub = [](Max left, Max right, Max carry_in = 0) -> bool {
        if(left < right)
            std::swap(left, right);

        Max highest_one = cast(Max) 1 << (BIT_SIZE<Max> - 1);
        Max adjusted_right = right & ~highest_one;
        Max adjusted_carry = cast(Max) low_bits(cast(T) carry_in);

        Max adjusted_left = left;
        if(left < adjusted_carry + adjusted_right)
            adjusted_left = add_no_overflow(left, adjusted_carry);

        Max expected = adjusted_left - adjusted_right - adjusted_carry;
        return test_sub_overflow_batch(adjusted_left, adjusted_right, adjusted_carry, Res{expected, 0});
    };

    //if constexpr(std::is_same_v<T, u8>)
    //{
    //    assert(test_sub_overflow_batch(0, 0, 0, Res{0, 0}));
    //    assert(test_sub_overflow_batch(1, 0, 0, Res{1, 0}));
    //    assert(test_sub_overflow_batch(0, 1, 0, Res{0xFF, 1}));
    //    assert(test_sub_overflow_batch(1, 1, 0, Res{0, 0}));
    //    assert(test_sub_overflow_batch(2, 1, 0, Res{1, 0}));

    //    assert(test_sub_overflow_batch(0xFF0000, 0x000001, 0, Res{0xFEFFFF, 0}));
    //    assert(test_sub_overflow_batch(0x000000, 0x000001, 0, Res{0xFF, 1}));
    //    assert(test_sub_overflow_batch(0xFF0000, 0x0000FF, 0, Res{0xFEFF01, 0}));

    //    assert(test_sub_overflow_batch(0xFFFF, 0xFF, 0, Res{0xFF00, 0}));
    //    assert(test_sub_overflow_batch(0xFFFE, 0xFF, 0, Res{0xFEFF, 0}));
    //    assert(test_sub_overflow_batch(0xFE, 0xFF, 0, Res{0xFF, 1}));
    //    assert(test_sub_overflow_batch(0x0, 0xFFFFFFFF, 0, Res{0x1, 1}));
    //    assert(test_sub_overflow_batch(0x1, 0xFFFFFFFF, 0, Res{0x2, 1}));
    //    assert(test_sub_overflow_batch(0xFFFFFFFF, 0x1, 0, Res{0xFFFFFFFE, 0}));
    //}

    //assert(auto_test_sub(0xFF, 0xFF));
    //assert(auto_test_sub(0xFFFF'FFFF'FFFF'FFFF, 0xFFFF'FFFF'FFFF'FFFF));
    //assert(auto_test_sub(56775421, 32101454));
    //assert(auto_test_sub(64344, 45));
    //assert(auto_test_sub(984643131, 81766));
    //assert(auto_test_sub(18748341324, 13875));
    //assert(auto_test_sub(23446634457799, 454406999));
    //assert(auto_test_sub(74984193, 982387980));
    assert(auto_test_sub(0x04b941000665a200, 0x7d57e, 0xc87afb6c));

    let distribution = make_max_distribution<T>();
    for(size_t i = 0; i < random_runs; i++)
    {
        Max left = distribution(*generator);
        Max right = distribution(*generator);
        Max carry = distribution(*generator);
        assert(auto_test_sub(left, right, carry));
    };
}

template <typename T>
runtime_proc test_complement_overflow(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    //we dont test oveflow of this op 
    let complement_overflow_batch = [](Max left_, Max carry_in = 1) -> Max{
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        mut pad_left = make_padded_big_int<T>(left_);

        Slice<T> out_s = out;
        Slice<T> left_s = pad_left.num;
        Slice<T> pad_left_s = pad_left.whole_num;

        let res1 = ::complement_overflow_batch<T>(&out_s, left_s, cast(T) carry_in);
        let res2 = ::complement_overflow_batch<T>(&pad_left_s, left_s, cast(T) carry_in);

        let num1 = unwrap(to_number(res1.slice));
        let num2 = unwrap(to_number(res2.slice));

        assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
        assert(num1 == num2 && "self assign must work correctly");

        return num1;
    };

    assert(complement_overflow_batch(0) == 0);
    assert(complement_overflow_batch(1) == cast(u8)(-1));
    assert(complement_overflow_batch(5) == cast(u8)(-5));
    assert(complement_overflow_batch(0xFF) == cast(u8)(-0xFF));
    assert(complement_overflow_batch(0xABCD) == cast(u16)(-0xABCD));
    assert(complement_overflow_batch(0xFFFF0000) == cast(u32)(-cast(i64)0xFFFF0000));
}

template <typename T>
runtime_proc test_shift_overflow(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    constexpr let shift_both = [](Max left_, Max right, Max carry_in, bool up_down, Iter_Direction direction) -> Batch_Op_Result{
        mut pad_left = make_padded_big_int<T>(left_, 1, 1);
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        Slice<T> out_s = out;
        Slice<T> left_s = pad_left.num;
        Slice<T> pad_left_s;

        if(direction == Iter_Direction::FORWARD)
            pad_left_s = pad_left.whole_num;
        else
            pad_left_s = slice(pad_left.whole_num, 2);

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

    constexpr let shift_up_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
        let forward = shift_both(left, right, carry_in, true, Iter_Direction::FORWARD);
        let backward = shift_both(left, right, carry_in, true, Iter_Direction::BACKWARD);

        assert(forward.value == backward.value);
        assert(forward.overflow == backward.overflow);

        return forward;
    };

    constexpr let shift_down_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
        let forward = shift_both(left, right, carry_in, false, Iter_Direction::FORWARD);
        let backward = shift_both(left, right, carry_in, false, Iter_Direction::BACKWARD);

        assert(forward.value == backward.value);
        assert(forward.overflow == backward.overflow);

        return forward;
    };

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
}

template <typename T>
runtime_proc test_mul_overflow(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    let mul_overflow_batch = [](Max left_, Max right, Max carry_in = 0) -> Batch_Op_Result{
        mut pad_left = make_padded_big_int<T>(left_);
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        
        Slice<T> out_s = out;
        Slice<T> left_s = pad_left.num;
        Slice<T> pad_left_s = pad_left.whole_num;
        let info = Optim_Info{};

        let res1 = ::mul_overflow_batch<T>(&out_s, left_s, cast(T) right, info, cast(T) carry_in);
        let res2 = ::mul_overflow_batch<T>(&pad_left_s, left_s, cast(T) right, info, cast(T) carry_in);

        let num1 = unwrap(to_number(res1.slice));
        let num2 = unwrap(to_number(res2.slice));

        assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
        assert(num1 == num2 && "self assign must work correctly");

        return Batch_Op_Result{num1, res1.overflow};
    };

    //@TODO: add automated tests
    assert((mul_overflow_batch(0, 0) == Res{0, 0}));
    assert((mul_overflow_batch(1, 0) == Res{0, 0}));
    assert((mul_overflow_batch(0, 1) == Res{0, 0}));
    assert((mul_overflow_batch(1, 1) == Res{1, 0}));
    assert((mul_overflow_batch(0, 1, 0xFF) == Res{0, 0xFF}));
    assert((mul_overflow_batch(25, 1) == Res{25, 0}));
    assert((mul_overflow_batch(25, 5) == Res{125, 0}));
    assert((mul_overflow_batch(0xFF, 0x10) == Res{0xF0, 0xF}));
    assert((mul_overflow_batch(0xA1, 0x11) == Res{0xB1, 0xA}));
    assert((mul_overflow_batch(0xABCDEF00, 1) == Res{0xABCDEF00, 0}));
    assert((mul_overflow_batch(0xABCDEF00, 0) == Res{0, 0}));
    assert((mul_overflow_batch(0xABCDEF00, 0xAB) == Res{0xC28EA500, 0x72}));
    assert((mul_overflow_batch(13745, 137) == Res{0xBBB9, 0x1C}));
}

template <typename T>
runtime_proc test_div_overflow_low(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using SRes = Overflow<T>;
    using Big = Big_Int_<T>;

    assert((div_overflow_low<T>(0x12, 0x0A, 0) == SRes{1, 0x08}));

    let div_overflow_low_batch = [](Max left_, Max right, Max carry_in = 0) -> Batch_Op_Result{
        constexpr size_t padding_before = 2;

        mut pad_left = make_padded_big_int<T>(left_, padding_before, MAX_TYPE_SIZE_FRACTION<T>);
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

        Slice<T> out_s = out;
        Slice<T> left_s = pad_left.num;
        Slice<T> pad_left_s = slice<T>(pad_left.whole_num, padding_before + 1);
        let info = Optim_Info{};

        assert(unwrap(to_number(left_s)) == left_);

        let res1 = ::div_overflow_low_batch<T>(&out_s, left_s, cast(T) right, info, cast(T) carry_in);
        let res2 = ::div_overflow_low_batch<T>(&pad_left_s, left_s, cast(T) right, info, cast(T) carry_in);

        let num1 = unwrap(to_number(res1.slice));
        let num2 = unwrap(to_number(res2.slice));

        assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
        assert(num1 == num2 && "self assign must work correctly");

        return Batch_Op_Result{num1, res1.overflow};
    };

    assert((div_overflow_low_batch(0, 1) == Res{0, 0}));
    assert((div_overflow_low_batch(1, 1) == Res{1, 0}));
    assert((div_overflow_low_batch(25, 1) == Res{25, 0}));
    assert((div_overflow_low_batch(25, 5) == Res{5, 0}));
    assert((div_overflow_low_batch(0xFF, 0x08) == Res{31, 7}));
    assert((div_overflow_low_batch(0xA1, 0x0F) == Res{10, 11}));
    assert((div_overflow_low_batch(0xABCDEF00, 1) == Res{0xABCDEF00, 0}));
    assert((div_overflow_low_batch(0xABCDEF00, 0x07) == Res{411771428, 4}));
    assert((div_overflow_low_batch(0xABCDEF00, 0x0F) == Res{192160000, 0}));
    assert((div_overflow_low_batch(13745, 0xB) == Res{1249, 6}));
}

template <integral T>
runtime_proc test_mul_quadratic(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big_Int = Big_Int_<T>;

    proc test_fused_mul_add = [](Max left_, Max right_, Max coef, Max expected, Optim_Info info = Optim_Info{}) -> bool{
        Big_Int left = left_;
        Big_Int right = right_;
        Big_Int output = make_big_int_with_size<T>(left.size);

        Slice<T> output_s = output;
        Slice<T> left_s = left;
        Slice<T> right_s = right;

        let res_out = ::fused_mul_add_overflow_batch<T>(&output_s, left_s, right_s, cast(T) coef, info, Op_Location::OUT_OF_PLACE);
        let res_in = ::fused_mul_add_overflow_batch<T>(&left_s, left_s, right_s, cast(T) coef, info, Op_Location::IN_PLACE);
        
        push(&output, res_out.overflow);
        push(&left, res_in.overflow);

        let num_out = unwrap(to_number<Max, T>(output));
        let num_in = unwrap(to_number<Max, T>(left));

        return expected == num_out && expected == num_in; 
    };

    proc auto_test_fused_mul_add = [](Max left, Max right, Max coef, Optim_Info info = Optim_Info{}) -> bool{
        return test_fused_mul_add(left, right, coef, left + coef*right, info);
    };

    proc test_mul_quadratic = [](Max left_, Max right_, Max expected, Optim_Info info = Optim_Info{}) -> bool{
        Big_Int left = left_;
        Big_Int right = right_;
        Big_Int res;
        Big_Int temp;

        resize(&temp, max(left.size, right.size) + 1);
        resize(&res, left.size + right.size + 1);

        Slice<T> res_s = res;
        Slice<T> temp_s = temp;
        Slice<const T> left_s = left;
        Slice<const T> right_s = right;

        ::mul_quadratic<T>(&res_s, &temp_s, left_s, right_s, info);
        let normal = unwrap(to_number(res_s));

        ::mul_quadratic_fused<T>(&res_s, left_s, right_s, info);
        let fused = unwrap(to_number(res_s));

        return expected == normal && expected == fused;
    };

    let auto_test_mul_quadratic = [=](Max left, Max right, Optim_Info info) -> bool
    {
        //to protect from overflow we check if the result can overflow by computing the necessary ammount of bits to represent the number
        // (by doing find_last_set_bit) and then shifting it so that they will not overflow

        size_t adjusted_left = left;
        size_t adjusted_right = right;

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
        
        return test_mul_quadratic(adjusted_left, adjusted_right, adjusted_left * adjusted_right, info);
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

    let distribution = make_max_distribution<T>();
    for(size_t i = 0; i < random_runs; i++)
    {
        Optim_Info info = generate_random_optims(generator);
        Max left = distribution(*generator);
        Max right = distribution(*generator);
        assert(auto_test_mul_quadratic(left, right, info));
    };
}

template <integral T>
runtime_proc test_div_bit_by_bit(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big_Int = Big_Int_<T>;

    enum Will_Fail
    {
        FAIL,
        OKAY
    };

    let test_div_bit_by_bit = [](Max left_, Max right_, Max ex_quotient, Max ex_remainder, Will_Fail expected_fail = OKAY, Optim_Info info = Optim_Info{}) -> bool{
        Big_Int left = left_; 
        Big_Int right = right_;
        Big_Int res;
        Big_Int rem;

        resize(&rem, max(left.size, right.size));
        resize(&res, max(left.size, right.size));

        Slice<T> res_s = res;
        Slice<T> rem_s = rem;
        Slice<const T> left_s = left;
        Slice<const T> right_s = right;

        let div_res = ::div_bit_by_bit<T>(&res_s, &rem_s, left_s, right_s, info);
        if(expected_fail == FAIL)
            return div_res.has == false;

        let quotient = unwrap(to_number(res_s));
        let remainder = unwrap(to_number(rem_s));
        return ex_quotient == quotient && ex_remainder == remainder && div_res.has;
    };


    let auto_test_div_bit_by_bit = [&](Max left, Max right, Optim_Info info = Optim_Info{}) -> bool{
        if(right == 0)
            return test_div_bit_by_bit(left, right, 0, 0, FAIL);

        Max quotient = left / right;
        Max remainder = left % right;
        return test_div_bit_by_bit(left, right, quotient, remainder, OKAY, info);
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

    let distribution = make_max_distribution<T>();
    for(size_t i = 0; i < random_runs; i++)
    {
        Optim_Info info = generate_random_optims(generator);
        Max_Unsigned_Type left = distribution(*generator);
        Max_Unsigned_Type right = distribution(*generator);
        assert(auto_test_div_bit_by_bit(left, right, info));
    };
}


template <typename T>
void println(Slice<T> slice)
{
    std::cout << "[";
    if(slice.size > 0)
    {
        std::cout << cast(Max_Unsigned_Type) slice[0];

        for(size_t i = 1; i < slice.size; i++)
            std::cout << ", " << cast(Max_Unsigned_Type) slice[i];
    }
    std::cout << "]\n";
}


template <integral From, integral To>
runtime_proc test_to_base(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    func test_to_base = [](Slice<const From> num, To base, Slice<const To> expected_rep, Optim_Info info) -> bool {
        Big_Int_<From> temp = Big_Int_<From>(num);
        Big_Int_<To> out = make_big_int_with_size<To>(required_size_to_base<From>(num.size, base));

        Slice<From> temp_s = temp;
        Slice<To> out_s = out;

        Slice<To> converted = (to_base<From, To>(&out_s,&temp_s, num, base, info));

        assert(is_reverse_striped_form(expected_rep));
        assert(is_reverse_striped_form(converted));

        return is_equal<To>(converted, expected_rep);
    };

    func manual_test_to_base = [](Max num, To base, std::initializer_list<To> expected_rep, Optim_Info info = Optim_Info{}) -> bool 
    {
        Big_Int_<From> num_ = Big_Int_<From>(num);
        Slice<const To> expected_rep_s = {std::data(expected_rep), std::size(expected_rep)};
        return test_to_base(num_, base, expected_rep_s, info);
    };

    func auto_to_base = [](Max num, To base, Optim_Info info) -> bool 
    {
        Big_Int_<From> num_ = num;
        Big_Int_<To> expected_rep = make_big_int_with_size<To>(required_size_to_base<From>(num_.size, base));
        size_t size = 0;
        for(Max curr_val = num; curr_val != 0; curr_val /= base, size ++)
        {
            Max rem = curr_val % base;
            expected_rep[size] = cast(To) rem;
        }

        Slice<To> expected_rep_s = {std::data(expected_rep), size};
        reverse(&expected_rep_s);
        return test_to_base(num_, base, expected_rep_s, info);
    };

    func test_from_base = [](Slice<const From> rep, From base, Slice<const To> expected_num, Optim_Info info) -> bool {
        size_t required_size = required_size_from_base<To>(rep.size, base);
        Big_Int_<To> out = make_big_int_with_size<To>(required_size);

        Slice<To> out_s = out;
        size_t from_size = sizeof(From);
        size_t to_size = sizeof(To);

        Slice<To> converted = (from_base<From, To>(&out_s, rep, base, info));

        
        assert(is_reverse_striped_form(expected_num));
        assert(is_reverse_striped_form(converted));

        //println(expected_num);
        //println(converted);

        return is_equal<To>(converted, expected_num);
    };

    func manual_test_from_base = [](std::initializer_list<From> rep, From base, Max expected_num, Optim_Info info = Optim_Info{}) -> bool 
    {
        Big_Int_<To> expected_num_ = expected_num;
        Slice<const From> rep_s = {std::data(rep), std::size(rep)};
        return test_from_base(rep_s, base, expected_num_, info);
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

    //assert((manual_test_from_base({1}, 2, 1)));
    //assert((manual_test_from_base({1,0}, 2, 2)));
    //assert((manual_test_from_base({1,0,0}, 2, 4)));
    //assert((manual_test_from_base({1,0,1}, 2, 5)));
    //assert((manual_test_from_base({1,0,1,0}, 2, 10)));

    //assert((manual_test_from_base({}, 9, 0)));
    //assert((manual_test_from_base({1}, 9, 1)));
    //assert((manual_test_from_base({2}, 9, 2)));
    assert((manual_test_from_base({1,0,8,7,0,5}, 9, 65453)));
    assert((manual_test_from_base({1,5,2,6,2,1,1}, 9, 844354)));
    assert((manual_test_from_base({8,4,4,3,5,4}, 10, 844354)));

    let num_dist = make_max_distribution<From>();
    constexpr size_t max_1 = std::numeric_limits<To>::max() >> HALF_BIT_SIZE<To>;
    constexpr size_t max_2 = std::numeric_limits<From>::max() >> HALF_BIT_SIZE<From>;

    let base_dist = make_max_distribution<From, To>(2, cast(To) min(max_1, max_2));
    for(size_t i = 0; i < random_runs; i++)
    {
        Optim_Info info = generate_random_optims(generator);
        Max_Unsigned_Type num = num_dist(*generator);
        To base = base_dist(*generator);
        assert(auto_to_base(num, base, info));
    };
}

runtime_proc run_typed_tests()
{
    test_misc();
    //test_add_overflow();
    //test_complement_overflow();
    //test_sub_overflow();
    //test_shift_overflow();
    //test_mul_overflow();
    //test_div_overflow_low();
}

template <integral T>
runtime_proc run_untyped_tests(Random_Generator* generator, size_t random_runs)
{
    test_complement_overflow<T>(generator, random_runs);
    test_sub_overflow<T>(generator, random_runs);
    test_shift_overflow<T>(generator, random_runs);
    test_mul_overflow<T>(generator, random_runs);
    test_div_overflow_low<T>(generator, random_runs);
    test_mul_quadratic<T>(generator, random_runs);
    test_div_bit_by_bit<T>(generator, random_runs);

    const size_t quarter_runs = random_runs / 4;
    test_to_base<T, u8>(generator, quarter_runs);
    test_to_base<T, u16>(generator, quarter_runs);
    test_to_base<T, u32>(generator, quarter_runs);
    test_to_base<T, u64>(generator, quarter_runs);
}

runtime_proc run_tests()
{
    std::random_device os_seed;
    const u32 seed = os_seed();

    Random_Generator generator(seed);
    const size_t random_runs = 10000;

    run_typed_tests();
    run_untyped_tests<u8>(&generator, random_runs);
    run_untyped_tests<u16>(&generator, random_runs);
    run_untyped_tests<u32>(&generator, random_runs);
    run_untyped_tests<u64>(&generator, random_runs);
}
