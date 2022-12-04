#pragma once


#include <iostream>
#include <cctype>
#include <random>
#include <limits>
#include <initializer_list>
#include <algorithm>
#include <array>

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
func make_sized_big_int(size_t out_size) -> Big_Int_<T>
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


template <integral T, integral Of = Max_Unsigned_Type>
runtime_func make_exponential_distribution(Of min = 0, Of max = std::numeric_limits<Of>::max()) -> Uniform_Exponential_Distribution<Of>
{
    double quot = pow(2, BIT_SIZE<T>);
    return {cast(double) min, cast(double) max, quot};
}

template <integral Of = Max_Unsigned_Type>
runtime_func make_uniform_distribution(Of min = 0, Of max = std::numeric_limits<Of>::max()) -> std::uniform_int_distribution<Of>
{
    return std::uniform_int_distribution<Of>{min, cast(double) max};
}

runtime_func generate_random_optims(Random_Generator* gen) -> Optim_Info
{
    std::uniform_int_distribution<> dist{0, 1};
    let random_bool = [&](){ return cast(bool) dist(*gen); };

    Optim_Info optims;
    optims.mul_shift = random_bool();
    optims.mul_consts = random_bool();
    optims.mul_half_bits = random_bool();

    optims.div_shift = random_bool();
    optims.div_consts = random_bool();
    optims.rem_optims = random_bool();
    return optims;
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


u64 ipow_native(u64 x, u64 n)
{
    if (x <= 1) 
        return x;

    u64 y = 1;
    for (; n != 0; n >>= 1, x *= x)
        if (n & 1)
            y *= x;

    return y;
}


u64 ipow(u64 x, u64 n)
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

u64 trivial_ipow(u64 x, u64 n)
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
        u64 power = ipow(r | 1, n); //slow
        u64 bit = power <= (x >> s);
        r |= bit;
    }

    return r;
}
u64 max_root_iters = 0;
u64 max_root_value = 0;
u64 max_root_power = 0;

u64 max_log_iters = 0;
u64 max_log_value = 0;
u64 max_log_power = 0;


u64 max_blog_iters = 0;
u64 max_blog_value = 0;
u64 max_blog_power = 0;

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

    u64 iters = 0;
    while(true)
    {
        // the newton update
        u64 new_upper_r = (n-1) * r + x / ipow(r, n-1);

        prev_r = r;
        r = new_upper_r / n;

        //we monotonically decrease 
        // => prev_r should always be bigger than r
        // => if not we break and return the last value that 
        //     satisfied the condition (prev_r)
        if(r >= prev_r)
            break;
        iters ++;
    }

    if(iters > max_root_iters)
    {
        max_root_iters = iters;
        max_root_value = of_value;
        max_root_power = power;
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
    u64 iters = 0;

    while(true)
    {
        u64 value_from_curr_approximate = ipow(b, y);
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
        iters ++;
    }


    if(iters > max_log_iters)
    {
        max_log_iters = iters;
        max_log_value = of_value;
        max_log_power = base;
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
    u64 iters = 0;

    while(true)
    {
        if(curr_hi - curr_lo <= 1)
            break;

        u64 value_from_curr_approximate = ipow(b, y);
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
        iters ++;
    }

    if(iters > max_blog_iters)
    {
        max_blog_iters = iters;
        max_blog_value = of_value;
        max_blog_power = base;
    }

    return curr_lo;
}

runtime_proc test_misc(Random_Generator* generator, size_t random_runs)
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


    assert(ipow(3,3) == 27);

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
        u64 powed_low = ipow(rooted, power);
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
            u64 ref_powed_low = trivial_ipow(rooted, power);
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

    std::cout << "ROOT" << std::endl;
    std::cout << "iters: " << max_root_iters << std::endl;
    std::cout << "value: " << max_root_value << std::endl;
    std::cout << "power: " << max_root_power << std::endl;

    std::cout << "LOG" << std::endl;
    std::cout << "iters: " << max_log_iters << std::endl;
    std::cout << "value: " << max_log_value << std::endl;
    std::cout << "power: " << max_log_power << std::endl;


    std::cout << "BLOG" << std::endl;
    std::cout << "iters: " << max_blog_iters << std::endl;
    std::cout << "value: " << max_blog_value << std::endl;
    std::cout << "power: " << max_blog_power << std::endl;
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
        mut out = make_sized_big_int<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

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

    proc test_single_add = [](Max left, Max single){
        Max highest_one = cast(Max) 1 << (BIT_SIZE<Max> - 1);
        Max adjusted_left = left & ~highest_one;
        T adjusted_single = cast(T) single & ~highest_one;
        Max expected = adjusted_left + adjusted_single;

        Big_Int_<T> left_copy = adjusted_left;
        Big_Int_<T> left_ = adjusted_left; 

        Slice<T> left_copy_s = left_copy;
        let res2 = add_overflow_batch<T>(&left_copy_s, left_copy_s, adjusted_single, Op_Location::IN_PLACE);

        bool ret = true;
        if(left != 0)
        {
            Slice<T> left_s = left_;
            let res1 = add_overflow_batch<T>(&left_s, left_s, Slice<T>{&adjusted_single, 1}, Op_Location::IN_PLACE);

            bool slices = is_equal<T>(res1.slice, res2.slice);
            bool overflows = res1.overflow == res2.overflow;

            ret = ret && slices && overflows;
        }
        else
        {
            assert(res2.overflow == adjusted_single);
        }

        return ret;
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
        Max left = distribution(*generator);
        Max right = distribution(*generator);
        Max carry = distribution(*generator);
        assert(auto_test_add(left, right, carry));
        assert(test_single_add(left, cast(T) right));
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
        mut out = make_sized_big_int<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

        Slice<T> out_s = out;
        Slice<T> pad_left_s = pad_left.whole_num;
        let loc = Op_Location::OUT_OF_PLACE;

        let res1 = ::sub_overflow_batch<T>(&out_s, left, right, loc, cast(T) carry_in);
        let res2 = ::sub_overflow_batch<T>(&pad_left_s, pad_left.num, right, loc, cast(T) carry_in); //test self assign

        Slice<T> res1_s = res1.slice;
        Slice<T> res2_s = res2.slice;

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

    proc auto_test_sub = [](Max left, Max right, Max carry_in = 0) -> bool {
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
    proc complement_overflow_batch = [](Max left_, Max carry_in = 1) -> Max{
        mut out = make_sized_big_int<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
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


    proc auto_test_complement = [](Max left, Max carry_in = 1) -> bool{
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
        Max left = exponential(*generator);
        Max carry = uniform(*generator);
        assert(auto_test_complement(left, carry));
    };
}

template <typename T>
runtime_proc test_shift_overflow(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    constexpr let shift_both = [](Max left_, Max right, Max carry_in, bool up_down, Iter_Direction direction) -> Batch_Op_Result{
        mut pad_left = make_padded_big_int<T>(left_, 1, 1);
        mut out = make_sized_big_int<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
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
}

struct Adjusted
{
    Max_Unsigned_Type left;
    Max_Unsigned_Type right;
};

runtime_func adjust_to_not_mul_overflow(Max_Unsigned_Type left, Max_Unsigned_Type right) -> Adjusted
{
    using Max = Max_Unsigned_Type;
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
runtime_proc test_mul_overflow(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big = Big_Int_<T>;

    proc test_mul_overflow_batch = [](Max left_, Max right, Max carry_in, Res expected, Optim_Info optims, bool check_overflow = true) -> bool {
        mut pad_left = make_padded_big_int<T>(left_);
        mut out = make_sized_big_int<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);
        
        Slice<T> out_s = out;
        Slice<T> left_s = pad_left.num;
        Slice<T> pad_left_s = pad_left.whole_num;

        let res1 = ::mul_overflow_batch<T>(&out_s, left_s, cast(T) right, optims, cast(T) carry_in);
        let res2 = ::mul_overflow_batch<T>(&pad_left_s, left_s, cast(T) right, optims, cast(T) carry_in);

        Slice<T> res1_s = res1.slice;
        Slice<T> res2_s = res2.slice;

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

    proc auto_test_mul = [](Max left, Max right, Optim_Info optims) -> bool
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
        }
    }

    assert(auto_test_mul(0x14156f, 0x3c, Optim_Info{}));

    let exponential = make_exponential_distribution<T, Max>();
    let uniform = std::uniform_int_distribution<long long>(0, FULL_MASK<T> >> 1);
    for(size_t i = 0; i < random_runs; i++)
    {
        Max left = exponential(*generator);
        Max right = uniform(*generator);
        let optims = generate_random_optims(generator);
        assert(auto_test_mul(left, right, optims));
    };
}

template <typename T>
runtime_proc test_div_overflow_low(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using SRes = Overflow<T>;
    using Big = Big_Int_<T>;

    assert((div_overflow_low<T>(0x12, 0x0A, 0) == SRes{1, 0x08}));

    proc test_div_overflow_low_batch = [](Max left_, Max right, Max carry_in, Res expected, Optim_Info optims) -> bool {
        constexpr size_t padding_before = 2;

        mut pad_left = make_padded_big_int<T>(left_, padding_before, MAX_TYPE_SIZE_FRACTION<T>);
        mut out = make_sized_big_int<T>(MAX_TYPE_SIZE_FRACTION<T> + 1);

        Slice<T> out_s = out;
        Slice<T> left_s = pad_left.num;
        Slice<T> pad_left_s = slice<T>(pad_left.whole_num, padding_before + 1);

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
    
    proc auto_test_div = [](Max left, Max right, Optim_Info optims) -> bool {
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
        }
    }

    let exponential = make_exponential_distribution<T, Max>();
    let uniform = std::uniform_int_distribution<long long>(0, low_mask<T>());
    for(size_t i = 0; i < random_runs; i++)
    {
        Max left = exponential(*generator);
        Max right = uniform(*generator);
        let optims = generate_random_optims(generator);
        assert(auto_test_div(left, right, optims));
    };
}

template <integral T>
runtime_proc test_mul_quadratic(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Res = Batch_Op_Result;
    using Big_Int = Big_Int_<T>;

    proc test_fused_mul_add = [](Max left_, Max right_, Max coef, Max expected, Optim_Info optims = Optim_Info{}) -> bool{
        Big_Int left = left_;
        Big_Int right = right_;
        Big_Int output = make_sized_big_int<T>(left.size);

        Slice<T> output_s = output;
        Slice<T> left_s = left;
        Slice<T> right_s = right;

        let res_out = ::fused_mul_add_overflow_batch<T>(&output_s, left_s, right_s, cast(T) coef, optims, Op_Location::OUT_OF_PLACE);
        let res_in = ::fused_mul_add_overflow_batch<T>(&left_s, left_s, right_s, cast(T) coef, optims, Op_Location::IN_PLACE);
        
        push(&output, res_out.overflow);
        push(&left, res_in.overflow);

        let num_out = unwrap(to_number<Max, T>(output));
        let num_in = unwrap(to_number<Max, T>(left));

        return expected == num_out && expected == num_in; 
    };

    proc auto_test_fused_mul_add = [](Max left, Max right, Max coef, Optim_Info optims = Optim_Info{}) -> bool{
        return test_fused_mul_add(left, right, coef, left + coef*right, optims);
    };

    proc test_mul_quadratic = [](Max left_, Max right_, Max expected, Optim_Info optims = Optim_Info{}) -> bool{
        Big_Int left = left_;
        Big_Int right = right_;
        Big_Int res;
        Big_Int aux;

        size_t to_size = required_mul_to_size(left.size, right.size);
        size_t aux_size1 = required_mul_quadratic_auxiliary_size(left.size, right.size);
        size_t aux_size2 = required_mul_karatsuba_auxiliary_size(left.size, right.size);

        size_t aux_size = max(aux_size1, aux_size2)*10; //to allow all karatsuba recursions

        resize(&aux, aux_size);
        resize(&res, to_size);

        Slice<T> res_s = res;
        Slice<T> aux_s = aux;
        Slice<const T> left_s = left;
        Slice<const T> right_s = right;

        ::mul_quadratic<T>(&res_s, &aux_s, left_s, right_s, optims);
        let normal = unwrap(to_number(res_s));

        ::mul_quadratic_fused<T>(&res_s, left_s, right_s, optims);
        let fused = unwrap(to_number(res_s));

        ::mul_karatsuba<T>(&res_s, &aux_s, left_s, right_s, optims);
        let karatsuba = unwrap(to_number(res_s));

        return expected == normal && expected == fused && expected == karatsuba;
    };

    proc auto_test_mul_quadratic = [](Max left, Max right, Optim_Info optims) -> bool
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
        Max left = distribution(*generator);
        Max right = distribution(*generator);

        assert(auto_test_mul_quadratic(left, right, optims));
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

    let test_div_bit_by_bit = [](Max left_, Max right_, Max ex_quotient, Max ex_remainder, Will_Fail expected_fail = OKAY, Optim_Info optims = Optim_Info{}) -> bool{
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

        let div_res = ::div_bit_by_bit<T>(&res_s, &rem_s, left_s, right_s, optims);
        if(expected_fail == FAIL)
            return div_res.has == false;

        let quotient = unwrap(to_number(res_s));
        let remainder = unwrap(to_number(rem_s));
        return ex_quotient == quotient && ex_remainder == remainder && div_res.has;
    };


    let auto_test_div_bit_by_bit = [&](Max left, Max right, Optim_Info optims = Optim_Info{}) -> bool{
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
        Max_Unsigned_Type left = distribution(*generator);
        Max_Unsigned_Type right = distribution(*generator);
        assert(auto_test_div_bit_by_bit(left, right, optims));
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

template <integral Num, integral Rep, typename Conv_Fn>
func to_base(Slice<const Num> num, Num base, Conv_Fn conv, Optim_Info optims) -> Big_Int_<Rep> {
    Big_Int_<Num> temp = Big_Int_<Num>(num);
    Big_Int_<Rep> out = make_sized_big_int<Rep>(required_size_to_base<Num>(num.size, base));

    Slice<Num> temp_s = temp;
    Slice<Rep> out_s = out;

    Slice<Rep> converted = (to_base<Num, Rep>(&out_s, &temp_s, num, base, conv, optims));
    assert(is_striped_representation(converted));

    resize(&out, converted.size);
    return out;
};

template <integral Num, integral Rep, typename Conv_Fn>
func from_base(Slice<const Rep> rep, Num base, Conv_Fn conv, Optim_Info optims) -> Big_Int_<Num> {
    size_t required_size = required_size_from_base<Num>(rep.size, base);
    Big_Int_<Num> out = make_sized_big_int<Num>(required_size + 1);

    Slice<Num> out_s = out;
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

template <integral Num = Max_Unsigned_Type>
func to_base(Slice<const Num> num, Num base = 10) -> Big_Int_<char> {
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

    return to_base<Num, char>(num, base, conversion, Optim_Info{});
}

template <integral Num = Max_Unsigned_Type>
func from_base(const char* str, Num base = 10) -> Trivial_Maybe<Big_Int_<Num>> {
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
    Big_Int_<Num> res = from_base<Num, char>(rep, base, conversion, Optim_Info{});

    if(ok)
        return wrap(res);
    else
        return {};
}

template <integral Num, integral Rep>
func native_to_base(Max_Unsigned_Type num, Num base) -> Big_Int_<Rep> 
{
    using Max = Max_Unsigned_Type;
    const size_t num_digits = digits_to_represent<Num, Max>();
    const size_t max_size = required_size_to_base<Num>(num_digits, base);

    Big_Int_<Rep> converted = make_sized_big_int<Rep>(max_size);
    size_t size = 0;
    for(Max curr_val = num; curr_val != 0; curr_val /= base, size ++)
    {
        Max rem = curr_val % base;
        converted[size] = cast(Rep) rem;
    }

    resize(&converted, size);
    Slice<Rep> converted_s = converted;
    reverse(&converted_s);

    return converted;
};

template <integral Num, integral Rep>
func native_from_base(Slice<const Rep> rep, Num base) -> Max_Unsigned_Type 
{
    using Max = Max_Unsigned_Type;
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

template <integral Num, integral Rep>
runtime_proc test_to_base(Random_Generator* generator, size_t random_runs)
{
    static_assert(sizeof(Num) >= sizeof(Rep));
    using Max = Max_Unsigned_Type;

    proc id_to_conversion = [](Num num) -> Rep {return cast(Rep) num;};
    proc id_from_conversion = [](Rep num) -> Num {return cast(Num) num;};

    func test_to_base = [](Slice<const Num> num, Num base, Slice<const Rep> expected_rep, Optim_Info optims) -> bool 
    {
        Big_Int_<Rep> converted = to_base<Num, Rep>(num, base, id_to_conversion, optims);
        
        assert(is_striped_representation(expected_rep));

        return is_equal<Rep>(converted, expected_rep);
    };

    func manual_test_to_base = [](Max num, Num base, std::initializer_list<Rep> expected_rep, Optim_Info optims = Optim_Info{}) -> bool 
    {
        Big_Int_<Num> num_ = Big_Int_<Num>(num);
        Slice<const Rep> expected_rep_s = {std::data(expected_rep), std::size(expected_rep)};
        return test_to_base(num_, base, expected_rep_s, optims);
    };

    func auto_test_to_base = [](Max num, Num base, Optim_Info optims) -> bool 
    {
        Big_Int_<Num> num_ = num;
        Big_Int_<Rep> expected_rep = native_to_base<Num, Rep>(num, base);
        return test_to_base(num_, base, expected_rep, optims);
    };

    func test_from_base = [](Slice<const Rep> rep, Num base, Slice<const Num> expected_num, Optim_Info optims) -> bool 
    {
        Big_Int_<Num> num = from_base<Num, Rep>(rep, base, id_from_conversion, optims);
        
        assert(is_striped_representation<const Num>(expected_num));
        assert(is_striped_representation<const Num>(num));

        return is_equal<Num>(num, expected_num);
    };

    func manual_test_from_base = [](std::initializer_list<Rep> rep, Num base, Max expected_num, Optim_Info optims = Optim_Info{}) -> bool 
    {
        Big_Int_<Num> expected_num_ = expected_num;
        Slice<const Rep> rep_s = {std::data(rep), std::size(rep)};
        return test_from_base(rep_s, base, expected_num_, optims);
    };

    func auto_test_from_base = [](Slice<const Rep> rep, Num base, Optim_Info optims) -> bool 
    {
        Max num = native_from_base<Num, Rep>(rep, base);
        Big_Int_<Num> expected_num = num;
        return test_from_base(rep, base, expected_num, optims);
    };

    func test_to_and_fro = [](Max value, Num base, Optim_Info optims) -> bool
    {
        Big_Int_<Num> initial = value;
        let rep = to_base<Num, Rep>(initial, base, id_to_conversion, optims);
        let num = from_base<Num, Rep>(rep, base, id_from_conversion, optims);
        Slice<const Num> num_s = num;
        let obtained_value = unwrap(to_number(num_s));

        bool match = is_equal<Num>(num, initial);
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
    
    for(size_t i = 0; i < random_runs; i++)
    {
        Optim_Info optims = generate_random_optims(generator);
        Max_Unsigned_Type num = num_dist(*generator);
        Num base = base_dist(*generator);

        assert(auto_test_to_base(num, base, optims));
        assert(test_to_and_fro(num, base, optims));
    };
}

template <integral T>
runtime_proc test_pow_by_squaring(Random_Generator* generator, size_t random_runs)
{
    using Max = Max_Unsigned_Type;
    using Big_Int = Big_Int_<T>;
    using CSlice = Slice<const T>;
    using Slice = Slice<T>;

    proc pow_by_squaring_ = [](CSlice num, T pow, Optim_Info optims) -> Big_Int {
        size_t required_size = required_pow_to_size(num, pow);
        size_t aux_size = required_pow_auxiliary_size(num, pow)*100; //for karatsuba
        Big_Int out = make_sized_big_int<T>(required_size);
        Big_Int aux = make_sized_big_int<T>(aux_size);

        Slice out_s = out;
        Slice aux_s = aux;

        Slice powed = ::pow_by_squaring<T>(&out_s, &aux_s, num, pow, optims);
        resize(&out, powed.size);

        return out;
    };

    proc pow_by_squaring = [](Max num, Max pow, Optim_Info optims = Optim_Info{}) -> Max {
        Big_Int num_ = num;
        Big_Int powed = pow_by_squaring_(num_, cast(T) pow, optims);
        Slice powed_s = powed;

        Max res = unwrap(to_number(powed_s));
        return res;
    };

    proc test_pow = [](CSlice num, T pow, CSlice expected, Optim_Info optims) -> bool {
        Big_Int powed = pow_by_squaring_(num, pow, optims);
        Slice powed_s = powed;

        assert(is_striped_number(powed_s));
        assert(is_striped_number(expected));
        
        return is_equal<T>(powed_s, expected);
    };

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

runtime_proc run_typed_tests(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
{
    test_misc(generator, random_runs);
}

template <integral T>
runtime_proc run_untyped_tests(Random_Generator* generator, size_t random_runs, size_t controlled_runs)
{
    test_complement_overflow<T>(generator, random_runs);
    test_add_overflow<T>(generator, random_runs);
    test_sub_overflow<T>(generator, random_runs);
    test_shift_overflow<T>(generator, random_runs);
    test_mul_overflow<T>(generator, random_runs, controlled_runs);
    test_div_overflow_low<T>(generator, random_runs, controlled_runs);
    test_mul_quadratic<T>(generator, random_runs);
    test_div_bit_by_bit<T>(generator, random_runs);
    test_pow_by_squaring<T>(generator, random_runs);


    const size_t quarter_runs = random_runs / 4;
    if constexpr(sizeof(T) >= sizeof(u8))
        test_to_base<T, u8>(generator, quarter_runs);

    if constexpr(sizeof(T) >= sizeof(u16))
        test_to_base<T, u16>(generator, quarter_runs);

    if constexpr(sizeof(T) >= sizeof(u32))
        test_to_base<T, u32>(generator, quarter_runs);

    if constexpr(sizeof(T) >= sizeof(u64))
        test_to_base<T, u64>(generator, quarter_runs);
}

runtime_proc run_tests()
{
    std::random_device os_seed;
    const u32 seed = os_seed();

    Random_Generator generator(seed);
    const size_t random_runs = 1000;
    const size_t controlled_runs = 50;

    run_typed_tests(&generator, random_runs, controlled_runs);
    run_untyped_tests<u8>(&generator, random_runs, controlled_runs);
    run_untyped_tests<u16>(&generator, random_runs, controlled_runs);
    run_untyped_tests<u32>(&generator, random_runs, controlled_runs);
    run_untyped_tests<u64>(&generator, random_runs, controlled_runs);
}
