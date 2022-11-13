#pragma once


#include <iostream>
#include <cctype>
#include <random>
#include <limits>


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

template <integral T>
proc test_add_carry()
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;

    {
        using Big = Big_Int_<u8>;

        let add_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
            
            mut out = make_big_int_with_size<u8>(MAX_TYPE_SIZE_FRACTION<u8> + 1);
            span<u8> out_s = out;

            let res = ::add_overflow_batch<u8>(&out_s, Big{left}, Big{right}, cast(u8) carry_in);
            return Batch_Op_Result{as_number(out_s), res.overflow};
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

    let auto_test_add = [](Max left, Max right, Max carry_in = 0) -> bool {
        using Big = Big_Int_<T>;
        let numeric_add_res = add_overflow<Max>(left, right, carry_in);
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        span<T> out_s = out;

        let big_add_res = ::add_overflow_batch<T>(&out_s, Big{left}, Big{right}, cast(T) carry_in);

        Max adjusted_overflow = big_add_res.overflow;
        Max adjusted_value = as_number(out_s);
        if(big_add_res.size != MAX_TYPE_SIZE_FRACTION<T>)
        {
            adjusted_value += cast(Max) big_add_res.overflow << (big_add_res.size * BIT_SIZE<T>);
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
}

proc test_sub_carry()
{
    using T = u8;
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    let sub_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        span<T> out_s = out;

        let res = ::sub_overflow_batch<T>(&out_s, Big{left}, Big{right}, cast(T) carry_in);
        return Batch_Op_Result{as_number(out_s), res.overflow};
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


proc test_complement_carry()
{
    using T = u8;
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    //we dont test oveflow of this op 
    let complement_overflow_batch = [](Max left, Max carry_in = 1) -> Max{
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        span<T> out_s = out;

        let res = ::complement_overflow_batch<T>(&out_s, Big{left}, cast(T) carry_in);
        return as_number(out_s);
    };

    assert(complement_overflow_batch(0) == 0);
    assert(complement_overflow_batch(1) == cast(u8)(-1));
    assert(complement_overflow_batch(5) == cast(u8)(-5));
    assert(complement_overflow_batch(0xFF) == cast(u8)(-0xFF));
    assert(complement_overflow_batch(0xABCD) == cast(u16)(-0xABCD));
    assert(complement_overflow_batch(0xFFFF0000) == cast(u32)(-cast(i64)0xFFFF0000));
}

proc test_mul_carry()
{
    using T = u8;
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    let mul_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        span<T> out_s = out;

        let res = ::mul_overflow_batch<T>(&out_s, Big{left}, cast(T) right, cast(T) carry_in);
        let num = as_number(out_s);
        return Batch_Op_Result{num, res.overflow};
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



proc test_div_carry()
{
    using T = u8;
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using SRes = Overflow<T>;
    using Big = Big_Int_<T>;

    assert((div_overflow<T>(0x12, 0x0A, 0) == SRes{1, 0x08}));
    assert((div_overflow<T>(0x12, 0xEA, 0) == SRes{0, 0x12}));
    assert((div_overflow2<T>(0x0, 0xFF, 0x01) == SRes{1, 0x01}));

    let div_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        span<T> out_s = out;

        let big_left = Big{left};
        let res = ::div_overflow_batch<T>(&out_s, big_left, cast(T) right, cast(T) carry_in);
        let num = as_number(out_s);
        return Batch_Op_Result{num, res.overflow};
    };

    //@TODO: add automated tests
    //assert((div_overflow_batch(0, 0) == Res{0, 0})); //do not test div by 0 yet
    //assert((div_overflow_batch(1, 0) == Res{0, 0}));
    assert((div_overflow_batch(0, 1) == Res{0, 0}));
    assert((div_overflow_batch(1, 1) == Res{1, 0}));
    assert((div_overflow_batch(25, 1) == Res{25, 0}));
    assert((div_overflow_batch(25, 5) == Res{5, 0}));
    assert((div_overflow_batch(0xFF, 0x10) == Res{0xF, 0xF}));
    assert((div_overflow_batch(0xA1, 0x11) == Res{0x9, 0x8}));
    assert((div_overflow_batch(0xABCDEF00, 1) == Res{0xABCDEF00, 0}));
    assert((div_overflow_batch(0xABCDEF00, 0xAB) == Res{0x101344C, 0x3C}));
    assert((div_overflow_batch(0xABCDEF00, 0x1000) == Res{0xABCDE, 0xF00}));
    assert((div_overflow_batch(13745, 137) == Res{100, 45}));
}


proc test_shift_carry()
{
    using T = u8;
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    let shift_up_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        span<T> out_s = out;

        let res = ::shift_up_overflow_batch<T>(&out_s, Big{left}, right, cast(T) carry_in);
        return Batch_Op_Result{as_number(out_s), res.overflow};
    };

    let shift_down_overflow_batch = [](Max left, Max right, Max carry_in = 0) -> Batch_Op_Result{
        mut out = make_big_int_with_size<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        span<T> out_s = out;

        let res = ::shift_down_overflow_batch<T>(&out_s, Big{left}, right, cast(T) carry_in);
        return Batch_Op_Result{as_number(out_s), res.overflow};
    };

    //[0xFF, 0xEE, 0xDD] [0xAA, 0xBB, 0xCC] >> 4
    //  LOWER               UPPER
    //              shift_up_overflow_batch(LOWER, 4, 0xDD) 
    //shift_up_overflow_batch(LOWER, 4, 0x00) = middle_shifted_out

    //[0x0F, 0xFE, 0xED] [0xE0 | 0x0A, 0x
                            

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
    assert((shift_down_overflow_batch(0b0101'0011, 7) == Res{0b0, 0b10100110}));
    assert((shift_down_overflow_batch(0b0101'0011, 7, 0b1100'0111) == Res{0b1000'1110, 0b10100110}));
    assert((shift_down_overflow_batch(0xFFEEDD, 4) == Res{0x0FFEED, 0xD0}));
    assert((shift_down_overflow_batch(0xFFEEDD, 4, 0x0A) == Res{0xAFFEED, 0xD0}));
    assert((shift_down_overflow_batch(0b0101010101110001, 5) == Res{0b0000001010101011, 0b10001000}));
    assert((shift_down_overflow_batch(0b0101010101110001, 5, 0b1010'1111) == Res{0b0111101010101011, 0b10001000}));
}

template <integral T>
proc test_mul_quadratic()
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big_Int = Big_Int_<T>;

    let test_mul_overflow_batch = [](Big_Int left, Big_Int right, Max expected) -> bool{
        Big_Int res;
        Big_Int temp;

        resize(&temp, max(left.size, right.size) + 1);
        resize(&res, left.size + right.size + 1);

        span<T> res_s = res;
        span<T> temp_s = temp;
        span<const T> left_s = left;
        span<const T> right_s = right;

        ::mul_quadratic<T>(&res_s, &temp_s, left_s, right_s);
        let normal = as_number(res_s);

        ::mul_quadratic_fused<T>(&res_s, left_s, right_s);
        let fused = as_number(res_s);

        return expected == normal && expected == fused;
    };


    let auto_test_mul_overflow_batch = [](Big_Int left, Big_Int right) -> bool{
        
    };

    assert(test_mul_overflow_batch(0, 0, 0));
    assert(test_mul_overflow_batch(1, 0, 0));
    assert(test_mul_overflow_batch(0, 1, 0));
    assert(test_mul_overflow_batch(1, 1, 1));
    assert(test_mul_overflow_batch(25, 1, 25));
    assert(test_mul_overflow_batch(25, 5, 125));
    assert(test_mul_overflow_batch(0xFF, 0x10, 0xFF0));
    assert(test_mul_overflow_batch(0xA1, 0x11, 0xAB1));
    assert(test_mul_overflow_batch(0xABCDEF00, 1, 0xABCDEF00));
    assert(test_mul_overflow_batch(0xABCDEF00, 0, 0));
    assert(test_mul_overflow_batch(13745, 137, 0x1CBBB9));
    assert(test_mul_overflow_batch(0xFFFF, 0x10000, 0xFFFF0000));
    assert(test_mul_overflow_batch(0xABCD, 0x10001, 0xABCDABCD));
    assert(test_mul_overflow_batch(0xABCDEF124, 0x1254, 0xC4CDA61BA7D0));
    assert(test_mul_overflow_batch(0x1254, 0xABCDEF124, 0xC4CDA61BA7D0));
}

template <integral T>
runtime_proc run_typed_tests()
{
    test_add_carry<T>();
    test_complement_carry();
    test_sub_carry();
    test_mul_carry();
    test_div_carry();
    test_shift_carry();
    test_mul_quadratic<T>();
}
runtime_proc run_tests()
{
    run_typed_tests<u8>();
    run_typed_tests<u16>();
    run_typed_tests<u32>();
    run_typed_tests<u64>();
}



//int main( void )
//{
//
//    std::uniform_int_distribution<Max_Unsigned_Type> distribution(0, std::numeric_limits<Max_Unsigned_Type>::max());
//    for(size_t i = 0; i < random_runs; i++)
//    {
//        Max_Unsigned_Type left = distribution(*generator);
//        Max_Unsigned_Type right = distribution(*generator);
//
//        Max_Unsigned_Type res = left + right;
//        Max_Unsigned_Type carry = cast(Max_Unsigned_Type) (res < left || res < right);
//
//
//    };
//
//    std::random_device os_seed;
//    const u32 seed = os_seed();
//
//    Random_Generator generator(seed);
//}