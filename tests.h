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
    using T = u32;
    using Res_S = Overflow<T>;
    using Max = Max_Unsigned_Type;

    assert(integral_log2<Max>(0) == 0);
    assert(integral_log2<Max>(1) == 0);
    assert(integral_log2<Max>(2) == 1);
    assert(integral_log2<Max>(3) == 1);
    assert(integral_log2<Max>(4) == 2);
    assert(integral_log2<Max>(7) == 2);
    assert(integral_log2<Max>(8) == 3);
    assert(integral_log2<Max>(76) == 6);
    assert(integral_log2<Max>(513) == 9);
    assert(integral_log2<Max>(641) == 9);
    assert(integral_log2<Max>(1024) == 10);
    assert(integral_log2<u64>(0xf4f1f63db9b7e800) == 63);

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

proc test_add_overflow()
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using T = u8;
    using Big = Big_Int_<T>;

    let add_overflow_batch = [](Max left_, Max right_, Max carry_in = 0) -> Batch_Op_Result{
        Big left = left_;
        Big right = right_;

        mut pad_left = make_padded_big_int<T>(left_);
        mut pad_right = make_padded_big_int<T>(right_);

        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        Slice<T> out_s = out;
        Slice<T> pad_left_s = pad_left.whole_num;
        Slice<T> pad_right_s = pad_right.whole_num;

        let res1 = ::add_overflow_batch<T>(&out_s, left, right, cast(T) carry_in);
        let res2 = ::add_overflow_batch<T>(&pad_left_s, pad_left.num, right, cast(T) carry_in);
        let res3 = ::add_overflow_batch<T>(&pad_right_s, left, pad_right.num, cast(T) carry_in);

        let num1 = unwrap(to_number(res1.slice));
        let num2 = unwrap(to_number(res2.slice));
        let num3 = unwrap(to_number(res3.slice));

        assert(res1.slice.size == res2.slice.size 
            && res1.slice.size == res3.slice.size 
            && "self assign must work correctly");
        assert(num1 == num2 && num1 == num3 && "self assign must work correctly");
        assert(res1.overflow == res2.overflow && res1.overflow == res3.overflow);
        return Batch_Op_Result{num1, res1.overflow};
    };

    assert((add_overflow_batch(0, 0) == Res{0, 0}));
    assert((add_overflow_batch(0, 1) == Res{1, 0}));
    assert((add_overflow_batch(1, 1) == Res{2, 0}));

    assert((add_overflow_batch(0x0025, 0x00FF) == Res{0x24, 1}));
    assert((add_overflow_batch(0x00FF, 0x0025) == Res{0x24, 1}));
    assert((add_overflow_batch(0xFFFFFFFF, 0xFFFFFFFF) == Res{0xFFFFFFFE, 1}));

    assert((add_overflow_batch(0x0025, 0x01FF) == Res{0x224, 0}));
    assert((add_overflow_batch(0x01FF, 0x0025) == Res{0x224, 0}));
    assert((add_overflow_batch(0x01FF, 0) == Res{0x01FF, 0}));
    assert((add_overflow_batch(0x01FF, 1) == Res{0x0200, 0}));
    assert((add_overflow_batch(0x01FF, 1, 1) == Res{0x0201, 0}));
    assert((add_overflow_batch(1, 0x01FF) == Res{0x0200, 0}));

    assert((add_overflow_batch(0x1225, 0xEFEF) == Res{0x0214, 1}));
}

proc test_sub_overflow()
{
    using T = u8;
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    let sub_overflow_batch = [](Max left_, Max right_, Max carry_in = 0) -> Batch_Op_Result{
        Big left = left_;
        Big right = right_;

        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        mut pad_left = make_padded_big_int<T>(left_);
        Slice<T> out_s = out;
        Slice<T> pad_left_s = pad_left.whole_num;

        let res1 = ::sub_overflow_batch<u8>(&out_s, left, right, cast(u8) carry_in);
        let res2 = ::sub_overflow_batch<u8>(&pad_left_s, pad_left.num, right, cast(u8) carry_in);

        let num1 = unwrap(to_number(res1.slice));
        let num2 = unwrap(to_number(res2.slice));

        assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
        assert(num1 == num2 && "self assign must work correctly");
        return Batch_Op_Result{num1, res1.overflow};
    };

    assert((sub_overflow_batch(0, 0) == Res{0, 0}));
    assert((sub_overflow_batch(1, 0) == Res{1, 0}));
    assert((sub_overflow_batch(0, 1) == Res{0xFF, 1}));
    assert((sub_overflow_batch(1, 1) == Res{0, 0}));
    assert((sub_overflow_batch(2, 1) == Res{1, 0}));

    assert((sub_overflow_batch(0xFF0000, 0x000001) == Res{0xFEFFFF, 0}));
    assert((sub_overflow_batch(0x000000, 0x000001) == Res{0xFF, 1}));
    assert((sub_overflow_batch(0xFF0000, 0x0000FF) == Res{0xFEFF01, 0}));

    assert((sub_overflow_batch(0xFFFF, 0xFF) == Res{0xFF00, 0}));
    assert((sub_overflow_batch(0xFFFE, 0xFF) == Res{0xFEFF, 0}));
    assert((sub_overflow_batch(0xFE, 0xFF) == Res{0xFF, 1}));
    assert((sub_overflow_batch(0x0, 0xFFFFFFFF) == Res{0x1, 1}));
    assert((sub_overflow_batch(0x1, 0xFFFFFFFF) == Res{0x2, 1}));
    assert((sub_overflow_batch(0xFFFFFFFF, 0x1) == Res{0xFFFFFFFE, 0}));
}

proc test_complement_overflow()
{
    using T = u8;
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

proc test_shift_overflow()
{
    using T = u8;
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
            res1 = ::shift_up_overflow_batch<T>(&out_s, left_s, right, cast(T) carry_in, direction);
            res2 = ::shift_up_overflow_batch<T>(&pad_left_s, left_s, right, cast(T) carry_in, direction);
        }
        else
        {
            res1 = ::shift_down_overflow_batch<T>(&out_s, left_s, right, cast(T) carry_in, direction);
            res2 = ::shift_down_overflow_batch<T>(&pad_left_s, left_s, right, cast(T) carry_in, direction);
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
    assert((shift_up_overflow_batch(0, 1, 0xFF) == Res{0, 0x1})); //this is sort of a weird edge case I dont know how to best handle
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
    assert((shift_down_overflow_batch(0, 1, 0xFF) == Res{0, 0x80})); //edge case
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

proc test_mul_overflow()
{
    using T = u8;
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    let mul_overflow_batch = [](Max left_, Max right, Max carry_in = 0) -> Batch_Op_Result{
        mut pad_left = make_padded_big_int<T>(left_);
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        
        Slice<T> out_s = out;
        Slice<T> left_s = pad_left.num;
        Slice<T> pad_left_s = pad_left.whole_num;

        let res1 = ::mul_overflow_batch<T>(&out_s, left_s, cast(T) right, cast(T) carry_in);
        let res2 = ::mul_overflow_batch<T>(&pad_left_s, left_s, cast(T) right, cast(T) carry_in);

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
    assert((mul_overflow_batch(25, 1) == Res{25, 0}));
    assert((mul_overflow_batch(25, 5) == Res{125, 0}));
    assert((mul_overflow_batch(0xFF, 0x10) == Res{0xF0, 0xF}));
    assert((mul_overflow_batch(0xA1, 0x11) == Res{0xB1, 0xA}));
    assert((mul_overflow_batch(0xABCDEF00, 1) == Res{0xABCDEF00, 0}));
    assert((mul_overflow_batch(0xABCDEF00, 0) == Res{0, 0}));
    assert((mul_overflow_batch(0xABCDEF00, 0xAB) == Res{0xC28EA500, 0x72}));
    assert((mul_overflow_batch(13745, 137) == Res{0xBBB9, 0x1C}));
}

proc test_div_overflow_low()
{
    using T = u8;
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

        assert(unwrap(to_number(left_s)) == left_);

        let res1 = ::div_overflow_low_batch<T>(&out_s, left_s, cast(T) right, cast(T) carry_in);
        let res2 = ::div_overflow_low_batch<T>(&pad_left_s, left_s, cast(T) right, cast(T) carry_in);

        let num1 = unwrap(to_number(res1.slice));
        let num2 = unwrap(to_number(res2.slice));

        assert(res1.slice.size == res2.slice.size && "self assign must work correctly");
        assert(num1 == num2 && "self assign must work correctly");

        return Batch_Op_Result{num1, res1.overflow};
    };

    //@TODO:
    //assert((div_overflow_low_batch(0, 0) == Res{0, 0})); //do not test div by 0 yet
    //assert((div_overflow_low_batch(1, 0) == Res{0, 0}));
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
runtime_proc test_untyped_add_carry(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;

    let auto_test_add = [](Max left, Max right, Max carry_in = 0) -> bool {
        using Big = Big_Int_<T>;
        let numeric_add_res = add_overflow<Max>(left, right, carry_in);
        mut pad_left = make_padded_big_int<T>(left, 2, MAX_TYPE_SIZE_FRACTION<T>);
        Slice<T> left_s = pad_left.num;
        Slice<T> out_s = pad_left.whole_num;

        let big_add_res = ::add_overflow_batch<T>(&out_s, left_s, Big{right}, cast(T) carry_in);

        Max adjusted_overflow = big_add_res.overflow;
        Max adjusted_value = unwrap(to_number(big_add_res.slice));
        if(big_add_res.slice.size != MAX_TYPE_SIZE_FRACTION<T>)
        {
            adjusted_value += cast(Max) big_add_res.overflow << (big_add_res.slice.size * BIT_SIZE<T>);
            adjusted_overflow = 0;
        }

        return numeric_add_res.value == adjusted_value && numeric_add_res.overflow == adjusted_overflow;
    };

    assert(auto_test_add(0xFF, 0xFF));
    assert(auto_test_add(0xFFFF'FFFF'FFFF'FFFF, 0xFFFF'FFFF'FFFF'FFFF));
    assert(auto_test_add(56775421, 32101454));
    assert(auto_test_add(64344, 45));
    assert(auto_test_add(984643131, 81766));
    assert(auto_test_add(18748341324, 13875));
    assert(auto_test_add(23446634457799, 454406999));
    assert(auto_test_add(74984193, 982387980));

    std::uniform_real_distribution<> ditsr1;
    mut ditsr2 = ditsr1;
    let distribution = make_max_distribution<T>();
    for(size_t i = 0; i < random_runs; i++)
    {
        Max_Unsigned_Type left = distribution(*generator);
        Max_Unsigned_Type right = distribution(*generator);
        assert(auto_test_add(left, right));
    };
}

template <integral T>
runtime_proc test_mul_quadratic(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big_Int = Big_Int_<T>;

    let test_mul_quadratic = [](Max left_, Max right_, Max expected) -> bool{
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

        ::mul_quadratic<T>(&res_s, &temp_s, left_s, right_s);
        let normal = unwrap(to_number(res_s));

        ::mul_quadratic_fused<T>(&res_s, left_s, right_s);
        let fused = unwrap(to_number(res_s));

        return expected == normal && expected == fused;
    };


    let auto_test_mul_quadratic = [=](Max left, Max right) -> bool
    {
        //to protect from overflow we check if the result can overflow by computing the necessary ammount of bits to represent the number
        // (by doing integral_log2) and then shifting it so that they will not overflow

        size_t adjusted_left = left;
        size_t adjusted_right = right;

        size_t log_left = integral_log2(left) + 1; 
        size_t log_right = integral_log2(right) + 1; 
        //integral_log2 returns position of highest set bit but we want "under which power of two is the whole number" 
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

        assert(integral_log2(adjusted_left) + integral_log2(adjusted_right) + 2 <= BIT_SIZE<Max>);
        
        return test_mul_quadratic(adjusted_left, adjusted_right, adjusted_left * adjusted_right);
    };


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
    assert(test_mul_quadratic(0xABCDEF124, 0x1254, 0xC4CDA61BA7D0));
    assert(test_mul_quadratic(0x1254, 0xABCDEF124, 0xC4CDA61BA7D0));

    let distribution = make_max_distribution<T>();
    for(size_t i = 0; i < random_runs; i++)
    {
        Max left = distribution(*generator);
        Max right = distribution(*generator);
        assert(auto_test_mul_quadratic(left, right));
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

    let test_div_bit_by_bit = [](Max left_, Max right_, Max ex_quotient, Max ex_remainder, Will_Fail expected_fail = OKAY) -> bool{
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

        let div_res = ::div_bit_by_bit<T>(&res_s, &rem_s, left_s, right_s);
        if(expected_fail == FAIL)
            return div_res.has == false;

        let quotient = unwrap(to_number(res_s));
        let remainder = unwrap(to_number(rem_s));
        return ex_quotient == quotient && ex_remainder == remainder && div_res.has;
    };


    let auto_test_div_bit_by_bit = [&](Max left, Max right) -> bool{
        if(right == 0)
            return test_div_bit_by_bit(left, right, 0, 0, FAIL);

        Max quotient = left / right;
        Max remainder = left % right;
        return test_div_bit_by_bit(left, right, quotient, remainder, OKAY);
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
        Max_Unsigned_Type left = distribution(*generator);
        Max_Unsigned_Type right = distribution(*generator);
        assert(auto_test_div_bit_by_bit(left, right));
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
    func test_to_base = [](Big_Int_<From> in, To base, Slice<const To> expected) -> bool {
        Big_Int_<From> temp = in;

        Big_Int_<To> out = make_big_int_with_size<To>(required_size_to_base<From>(in.size, base));

        Slice<From> in_s = in;
        Slice<From> temp_s = temp;
        Slice<To> out_s = out;

        Slice<To> converted = (to_base<From, To>(in_s, &temp_s, &out_s, base));

        assert(is_striped_form(expected));
        assert(is_striped_form(converted));

        return compare<To>(converted, expected) == 0;
    };

    func manual_test_to_base = [](Max num, To base, std::initializer_list<To> expected) -> bool 
    {
        Big_Int_<From> in = num;
        Slice<const To> expected_s = {std::data(expected), std::size(expected)};
        return test_to_base(std::move(in), base, expected_s);
    };

    func auto_to_base = [](Max num, To base) -> bool 
    {
        Big_Int_<From> in = num;
        Big_Int_<To> expected = make_big_int_with_size<To>(required_size_to_base<From>(in.size, base));
        size_t size = 0;
        for(Max curr_val = num; curr_val != 0; curr_val /= base, size ++)
        {
            Max rem = curr_val % base;
            expected[size] = cast(To) rem;
        }

        Slice<const To> expected_s = {std::data(expected), size};
        assert(is_striped_form(expected_s));

        return test_to_base(std::move(in), base, expected_s);
    };

    func test_from_base = [](Big_Int_<From> rep, To base, Slice<const To> expected) -> bool {
        Big_Int_<To> out = make_big_int_with_size<To>(required_size_from_base<From>(rep.size, base));

        Slice<From> rep_s = rep;
        Slice<To> out_s = out;

        Slice<To> converted = (from_base<From, To>(rep_s, &out_s, base));

        assert(is_striped_form(expected));
        assert(is_striped_form(converted));

        return compare<To>(converted, expected) == 0;
    }

    func manual_test_from_base = [](std::initializer_list<To> from, To base, Max expected) -> bool 
    {
        Big_Int_<From> in = expected;
        Slice<const To> expected_s = {std::data(expected), std::size(expected)};
        return test_to_base(std::move(in), base, expected_s);
    };


    assert((manual_test_to_base(1, 2, {1})));
    assert((manual_test_to_base(2, 2, {0,1})));
    assert((manual_test_to_base(4, 2, {0,0,1})));
    assert((manual_test_to_base(5, 2, {1,0,1})));
    assert((manual_test_to_base(10, 2, {0,1,0,1})));

    assert((manual_test_to_base(0, 9, {})));
    assert((manual_test_to_base(1, 9, {1})));
    assert((manual_test_to_base(2, 9, {2})));
    assert((manual_test_to_base(65453, 9, {5,0,7,8,0,1})));
    assert((manual_test_to_base(844354, 9, {1,1,2,6,2,5,1})));

    assert((manual_test_from_base({1}, 1, 2)));
    assert((manual_test_from_base({0,1}, 2, 2)));
    assert((manual_test_from_base({0,0,1}, 2, 4)));
    assert((manual_test_from_base({1,0,1}, 2, 5)));
    assert((manual_test_from_base({0,1,0,1}, 2, 10)));

    assert((manual_test_from_base({}, 9, 0)));
    assert((manual_test_from_base({1}, 9, 1)));
    assert((manual_test_from_base({2}, 9, 2)));
    assert((manual_test_from_base({5,0,7,8,0,1}, 9, 65453)));
    assert((manual_test_from_base({1,1,2,6,2,5,1}, 9, 844354)));

    let num_dist = make_max_distribution<From>();
    constexpr let max_1 = std::numeric_limits<To>::max() >> HALF_BIT_SIZE<To>;
    constexpr let max_2 = std::numeric_limits<From>::max() >> HALF_BIT_SIZE<From>;

    let base_dist = make_max_distribution<From, To>(2, min(max_1, max_2));
    for(size_t i = 0; i < random_runs; i++)
    {
        Max_Unsigned_Type num = num_dist(*generator);
        To base = base_dist(*generator);
        assert(auto_to_base(num, base));
    };
}

runtime_proc run_typed_tests()
{
    test_misc();
    test_add_overflow();
    test_complement_overflow();
    test_sub_overflow();
    test_shift_overflow();
    test_mul_overflow();
    test_div_overflow_low();
}

template <integral T>
runtime_proc run_untyped_tests(Random_Generator* generator, size_t random_runs)
{
    test_untyped_add_carry<T>(generator, random_runs);
    test_mul_quadratic<T>(generator, random_runs);
    test_div_bit_by_bit<T>(generator, random_runs);
    test_to_base<T, u8>(generator, random_runs);
    test_to_base<T, u16>(generator, random_runs);
    test_to_base<T, u32>(generator, random_runs);
    test_to_base<T, u64>(generator, random_runs);
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
