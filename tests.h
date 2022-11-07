#pragma once


#include <iostream>
#include <cctype>
#include <random>
#include <limits>


#include "benchmark.h"
#include "big_int.h"
#include "big_int_ops.h"


using Random_Generator = std::mt19937;


proc test_carry_ops(Max_Unsigned_Type left, Max_Unsigned_Type right, Max_Unsigned_Type expected, Max_Unsigned_Type expected_carry, let carry_op, Max_Unsigned_Type base_carry)
{
    using T = u8;
    using Big_Int = Big_Int_<T>;
    size_t type_size_ratio = sizeof(Max_Unsigned_Type) / sizeof(T);
    Big_Int_<T> l = left;
    Big_Int_<T> r = right;
    Big_Int_<T> res;
    resize(&res, type_size_ratio);

    span<T> res_s = res;
    span<const T> l_s = l;
    span<const T> r_s = r;
    assert(as_number(l_s) == left);
    assert(as_number(r_s) == right);

    let carry = carry_op(&res_s, l_s, r_s, cast(T) base_carry);
    let obtained = as_number(res_s);

    return carry == expected_carry && obtained == expected;
}

proc test_add_carry(Max_Unsigned_Type left, Max_Unsigned_Type right, Max_Unsigned_Type expected, Max_Unsigned_Type expected_carry)
{
    return test_carry_ops(left, right, expected, expected_carry, add_carry<u8>, 0);
}


proc test_complement_carry(Max_Unsigned_Type left, Max_Unsigned_Type expected)
{
    using T = u8;
    using Big_Int = Big_Int_<T>;
    size_t type_size_ratio = sizeof(Max_Unsigned_Type) / sizeof(T);
    Big_Int l = left;
    Big_Int res;
    resize(&res, type_size_ratio);

    span<T> res_s = res;
    span<const T> l_s = l;

    let carry = twos_complement_carry(&res_s, l_s);
    let obtained = as_number(res_s);

    return obtained == expected;
}

proc test_sub_carry(Max_Unsigned_Type left, Max_Unsigned_Type right, Max_Unsigned_Type expected, Max_Unsigned_Type expected_carry)
{
    return test_carry_ops(left, right, expected, expected_carry, sub_carry<u8>, 0);
}

proc test_mul_carry(Max_Unsigned_Type left, Max_Unsigned_Type right, Max_Unsigned_Type expected, Max_Unsigned_Type expected_carry, size_t size)
{
    using T = u8;
    using Big_Int = Big_Int_<T>;
    size_t type_size_ratio = sizeof(Max_Unsigned_Type) / sizeof(T);
    Big_Int l = left;
    Big_Int out;
    resize(&out, type_size_ratio);

    span<T> res_s = out;

    let res = mul_carry<T>(&res_s, l, cast(T) right);
    let obtained = as_number(res_s);
    return res.carry == expected_carry && obtained == expected && res.size == size;
}

proc test_mul_quadratic(Max_Unsigned_Type left, Max_Unsigned_Type right, Max_Unsigned_Type expected)
{
    using T = u8;
    using Big_Int = Big_Int_<T>;
    Big_Int l = left;
    Big_Int r = right;
    Big_Int res;
    Big_Int temp;

    resize(&temp, max(l.size, r.size) + 1);
    resize(&res, l.size + r.size + 1);

    span<T> res_s = res;
    span<T> temp_s = temp;
    span<const T> l_s = l;
    span<const T> r_s = r;
    assert(as_number(l_s) == left);
    assert(as_number(r_s) == right);

    mul_quadratic(&res_s, &temp_s, l_s, r_s);
    let obtained = as_number(res_s);

    return obtained == expected;
}

proc test_add_carry()
{
    assert(test_add_carry(0, 0, 0, 0));
    assert(test_add_carry(0, 1, 1, 0));
    assert(test_add_carry(1, 1, 2, 0));

    assert(test_add_carry(0x0025, 0x00FF, 0x24, 1));
    assert(test_add_carry(0x00FF, 0x0025, 0x24, 1));
    assert(test_add_carry(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE, 1));

    assert(test_add_carry(0x0025, 0x01FF, 0x224, 0));
    assert(test_add_carry(0x01FF, 0x0025, 0x224, 0));
    assert(test_add_carry(0x01FF, 0, 0x01FF, 0));
    assert(test_add_carry(0x01FF, 1, 0x0200, 0));
    assert(test_add_carry(1, 0x01FF, 0x0200, 0));

    assert(test_add_carry(0x1225, 0xEFEF, 0x0214, 1));
}

proc test_complement_carry()
{
    assert(test_complement_carry(0, 0));
    assert(test_complement_carry(1, cast(u8)(-1)));
    assert(test_complement_carry(5, cast(u8)(-5)));
    assert(test_complement_carry(0xFF, cast(u8)(-0xFF)));
    assert(test_complement_carry(0xABCD, cast(u16)(-0xABCD)));
    assert(test_complement_carry(0xFFFF0000, cast(u32)(-cast(i64)0xFFFF0000)));
}

proc test_sub_carry()
{
    assert(test_sub_carry(0, 0, 0, 0));
    assert(test_sub_carry(1, 0, 1, 0));
    assert(test_sub_carry(0, 1, 0xFF, 1));
    assert(test_sub_carry(1, 1, 0, 0));
    assert(test_sub_carry(2, 1, 1, 0));

    assert(test_sub_carry(0xFF0000, 0x000001, 0xFEFFFF, 0));
    assert(test_sub_carry(0x000000, 0x000001, 0xFF, 1));
    assert(test_sub_carry(0xFF0000, 0x0000FF, 0xFEFF01, 0));

    
    assert(test_sub_carry(0xFFFF, 0xFF, 0xFF00, 0));
    assert(test_sub_carry(0xFFFE, 0xFF, 0xFEFF, 0));
    assert(test_sub_carry(0xFE, 0xFF, 0xFF, 1));
    assert(test_sub_carry(0x0, 0xFFFFFFFF, 0x1, 1));
    assert(test_sub_carry(0x1, 0xFFFFFFFF, 0x2, 1));
    assert(test_sub_carry(0xFFFFFFFF, 0x1, 0xFFFFFFFE, 0));
}

proc test_mul_carry()
{
    assert(test_mul_carry(0, 0, 0, 0, 0));
    assert(test_mul_carry(1, 0, 0, 0, 0));
    assert(test_mul_carry(0, 1, 0, 0, 0));
    assert(test_mul_carry(1, 1, 1, 0, 1));
    assert(test_mul_carry(25, 1, 25, 0, 1));
    assert(test_mul_carry(25, 5, 125, 0, 1));
    assert(test_mul_carry(0xFF, 0x10, 0xF0, 0xF, 1));
    assert(test_mul_carry(0xA1, 0x11, 0xB1, 0xA, 1));
    assert(test_mul_carry(0xABCDEF00, 1, 0xABCDEF00, 0, 4));
    assert(test_mul_carry(0xABCDEF00, 0, 0, 0, 0));
    assert(test_mul_carry(13745, 137, 0xBBB9, 0x1C, 2));
}

proc test_mul_quadratic()
{
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
}

runtime_proc run_tests()
{
    test_add_carry();
    test_complement_carry();
    test_sub_carry();
    test_mul_carry();
    test_mul_quadratic();
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