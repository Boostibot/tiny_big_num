#define proc auto
#define func [[nodiscard]] auto
#define cast(...) (__VA_ARGS__)

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <cstring>

//Defaults:
#ifndef TINY_NUM_MAX_UNSIGNED_TYPE
    #define TINY_NUM_MAX_UNSIGNED_TYPE uint64_t
#endif 

#ifndef TINY_NUM_TYPE_DEF
    #define TINY_NUM_TYPE_DEF uint64_t
#endif 

#ifdef INCLUDED_NEW_TYPE
    #undef INCLUDED_NEW_TYPE
    #undef INCLUDED_TINY_NUM_TYPE_DEF
#endif

//if TINY_NUM_TYPE_1 is first time included (TINY_NUM_TYPE_1_INCLUDED wasnt defined yet casuse is defined after)
// => define unlock so that the library code gets compiledw with the type
#if defined(TINY_NUM_TYPE_1) && !defined(TINY_NUM_TYPE_1_INCLUDED)
    #define INCLUDED_NEW_TYPE TINY_NUM_TYPE_1

#elif defined(TINY_NUM_TYPE_2) && !defined(TINY_NUM_TYPE_2_INCLUDED)
    #define INCLUDED_NEW_TYPE TINY_NUM_TYPE_2

#elif defined(TINY_NUM_TYPE_3) && !defined(TINY_NUM_TYPE_3_INCLUDED)
    #define INCLUDED_NEW_TYPE TINY_NUM_TYPE_3

#elif defined(TINY_NUM_TYPE_4) && !defined(TINY_NUM_TYPE_4_INCLUDED)
    #define INCLUDED_NEW_TYPE TINY_NUM_TYPE_4

#elif defined(TINY_NUM_TYPE_DEF) && !defined(TINY_NUM_TYPE_DEF_INCLUDED) && !defined(INCLUDED_NEW_TYPE)
    #define INCLUDED_NEW_TYPE TINY_NUM_TYPE_DEF
    #define INCLUDED_TINY_NUM_TYPE_DEF
#endif

//define as included all that were included so far
#ifdef TINY_NUM_TYPE_DEF
    #define TINY_NUM_TYPE_DEF_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_1
    #define TINY_NUM_TYPE_1_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_2
    #define TINY_NUM_TYPE_2_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_3
    #define TINY_NUM_TYPE_3_INCLUDED
#endif

#ifdef TINY_NUM_TYPE_4
    #define TINY_NUM_TYPE_4_INCLUDED
#endif

//Optim switches:
#ifndef DO_OPTIM_MUL_SHIFT
    #define DO_OPTIM_MUL_SHIFT false
#endif

#ifndef DO_OPTIM_MUL_CONSTS
    #define DO_OPTIM_MUL_CONSTS false
#endif

#ifndef DO_OPTIM_MUL_HALF_BITS
    #define DO_OPTIM_MUL_HALF_BITS false
#endif

#ifndef DO_OPTIM_DIV_SHIFT
    #define DO_OPTIM_DIV_SHIFT false
#endif

#ifndef DO_OPTIM_DIV_CONSTS
    #define DO_OPTIM_DIV_CONSTS false
#endif

#ifndef DO_OPTIM_REM_OPTIMS
    #define DO_OPTIM_REM_OPTIMS false
#endif

#ifndef MAX_RECURSION_DEPTH 
    #define MAX_RECURSION_DEPTH 20
#endif

#ifndef MUL_QUADRATIC_SINGLE_BELOW_SIZE 
    #define MUL_QUADRATIC_SINGLE_BELOW_SIZE 3
#endif

#ifndef MUL_QUADRATIC_BOTH_BELOW_SIZE 
    #define MUL_QUADRATIC_BOTH_BELOW_SIZE 4
#endif

#ifndef TRIVIAL_POW_BELOW_POWER 
    #define TRIVIAL_POW_BELOW_POWER 4
#endif



//Library instance wide helpers
#ifndef TINY_NUM_FIRST_TIME_INCLUDED
#define TINY_NUM_FIRST_TIME_INCLUDED
namespace helpers
{
    enum class State
    {
        OK = 0,
        MATH_ERROR,
        OUT_SIZE_TOO_SMALL,
        AUX_SIZE_TOO_SMALL,
        NUMBER_TOO_BIG, //some ops require number to be only low bits
        //NUMBER_NOT_STRIPPED, //considered hard errors that should be prevented! (requires outside manipulation to get into invalid state)
        //INVALID_ALIASING,
    };

    enum class Location
    {
        IN_PLACE,
        OUT_OF_PLACE
    };

    enum class Iter_Direction 
    {
        FORWARD,
        BACKWARD
    };

    struct Optims
    {
        bool mul_shift = DO_OPTIM_MUL_SHIFT;
        bool mul_consts = DO_OPTIM_MUL_CONSTS;
        bool mul_half_bits = DO_OPTIM_MUL_HALF_BITS;
        bool div_shift = DO_OPTIM_DIV_SHIFT;
        bool div_consts = DO_OPTIM_DIV_CONSTS;
        bool rem_optims = DO_OPTIM_REM_OPTIMS;

        size_t max_recursion_depth = MAX_RECURSION_DEPTH;
        size_t mul_quadratic_both_below_size = MUL_QUADRATIC_BOTH_BELOW_SIZE;
        size_t mul_quadratic_single_below_size = MUL_QUADRATIC_SINGLE_BELOW_SIZE;
        size_t pow_trivial_below_power = TRIVIAL_POW_BELOW_POWER;
    };  

    template <typename T>
    proc swap(T* a, T* b)
    {
        T temp = *a;
        *a = *b;
        *b = temp;
    }

    func max(size_t a, size_t b) -> size_t
    {
        return a > b ? a : b;
    }

    func min(size_t a, size_t b) -> size_t
    {
        return a < b ? a : b;
    }

    func div_round_up(size_t value, size_t to_multiple_of) -> size_t
    {
        return (value + to_multiple_of - 1) / to_multiple_of;
    }

    func pop_count_u32(uint32_t val) -> size_t
    {
        uint32_t i = cast(uint32_t) val;
        i = i - ((i >> 1) & 0x55555555);                    // add pairs of bits
        i = (i & 0x33333333) + ((i >> 2) & 0x33333333);     // quads
        i = (i + (i >> 4)) & 0x0F0F0F0F;                    // groups of 8
        return cast(size_t) (i * 0x01010101) >> 24;         // horizontal sum of bytes
    }

    using umax = uint64_t;

    func find_last_set_bit(umax val) -> size_t  
    {
        assert(sizeof(umax) <= 8 && "Higher integer sizes not supported");
        size_t  k = 0;
        if (val > 0xFFFFFFFFu) { val >>= 32; k |= 32; }
        if (val > 0x0000FFFFu) { val >>= 16; k |= 16; }
        if (val > 0x000000FFu) { val >>= 8;  k |= 8;  }
        if (val > 0x0000000Fu) { val >>= 4;  k |= 4;  }
        if (val > 0x00000003u) { val >>= 2;  k |= 2;  }

        k |= (val & 2) >> 1;
        return k;
    }

    func log2(umax val) -> size_t  
    {
        return find_last_set_bit(val);
    }

    func find_first_set_bit(umax val) -> size_t  
    {
        if(val == 0)
            return 0;

        //10110000
        //&   1111 -> nothing => +4 and shift
        //00001011
        //&     11 -> something => +0
        //00001011
        //&      1 -> somehing => +0
        // => answer: 4 + 0 + 0 = 4 

        assert(sizeof(umax) <= 8 && "Higher integer sizes not supported");
        size_t  k = 0;

        if ((val & 0xFFFFFFFFu) == 0) { val >>= 32; k |= 32; }
        if ((val & 0x0000FFFFu) == 0) { val >>= 16; k |= 16; }
        if ((val & 0x000000FFu) == 0) { val >>= 8;  k |= 8;  }

        if ((val & 0x0000000Fu) == 0) { val >>= 4;  k |= 4;  }
        if ((val & 0x00000003u) == 0) { val >>= 2;  k |= 2;  }

        k |= (~val & 0b1);
        return k;
    }

    func pop_count(umax val) -> size_t
    {
        uint32_t low = cast(uint32_t) (val & 0xFFFF'FFFF);
        uint32_t high = cast(uint32_t) (val >> 32);

        return pop_count_u32(low) + pop_count_u32(high);
    }

    func single_pow_trivial(umax x, umax n) -> umax
    {
        umax res = 1;
        for(umax i = 0; i < n; i++)
            res *= x;

        return res;
    }

    func single_pow_by_squaring(umax x, umax n) -> umax
    {
        if (x <= 1) 
            return x;
        if(n == 0)  
            return 1;

        umax i = 0;
        umax y = 1;
        size_t max_pos = find_last_set_bit(n);
        for(; i <= max_pos; i++)
        {
            umax bit = n & (1ull << i);
            if(bit)
                y *= x;
            x *= x;
        }
        return y;
    }

    func single_root_shifting(umax x, umax n) -> umax
    {
        if(n == 0)
            return 1;
        if (x <= 1) 
            return x;

        umax r = 1;
        umax s = cast(umax) ((find_last_set_bit(cast(umax) x) / n) * n);

        while(true)
        {
            if(s < n)
                break;

            s -= n; 
            r <<= 1; 
            umax power = single_pow_by_squaring(r | 1, n);
            umax bit = power <= (x >> s);
            r |= bit;
        }

        return r;
    }

    func single_root_newton(umax of_value, umax power) -> umax
    {
        const umax x = of_value; 
        const umax n = power;

        if(x == 0)
            return 1;
        if (x <= 1) 
            return x;

        //umax u = x;
        //umax s = x+1;

        // We want to find such r that:
        //   r^n = x
        // this can be rewritten as: 
        //   n*log(r) = log(x)
        // factoring r gives:
        //   r = e^(log(x) / n)
        //   r = 2^(log2(x) / n)

        // we want to find an initial estimate above such r so we need to increase the expression:
        //  we will perform this by parts:
        //  log2(x) <= [log2(x)] + 1
        //  log2(x)/n <= upper[ ([log2(x)] + 1) / n ] = [ ([log2(x)] + 1 + n - 1) / n ]
        //            = [ ([log2(x)] + n) / n ]
        // 
        //  so the upper estimate for r is:
        //      2^[ ([log2(x)] + n) / n ]

        const umax log2x = cast(umax) find_last_set_bit(cast(umax) x);

        const umax upper_initial_estimate = cast(umax) 1 << ((log2x + n) / n);
        umax r = upper_initial_estimate;
        umax prev_r = -1;

        while(true)
        {
            // the newton update
            umax new_upper_r = (n-1) * r + x / single_pow_by_squaring(r, n-1);

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

    func single_log_heyley(umax of_value, umax base) -> umax
    {
        const umax x = of_value;
        const umax b = base;

        //b^y = x
        //y = log_b(x)

        if(x <= 1)
            return 0;
        if(base <= 1)
            return 0;

        const umax log2x = cast(umax) find_last_set_bit(x);
        const umax log2b = cast(umax) find_last_set_bit(b);

        assert(log2x != 0);
        assert(log2b != 0);

        const umax lower_initial_estimate = log2x / (log2b + 1);
        const umax higher_initial_estimate = (log2x + log2b) / (log2b);
        umax y = higher_initial_estimate;
        umax prev_y = y;

        while(true)
        {
            umax value_from_curr_approximate = single_pow_by_squaring(b, y);
            umax c_x = value_from_curr_approximate;

            prev_y = y;
            if(c_x <= x)
                break;

            umax new_y = y - 2*(c_x - x) / (x + c_x);

            //if we didnt move at all we are one above the answer
            if(new_y == y)
                return new_y - 1;

            y = new_y;
        }

        return prev_y;

    }

    func single_log_bsearch(umax of_value, umax base) -> umax
    {
        const umax x = of_value;
        const umax b = base;

        if(x <= 1)
            return 0;
        if(base <= 1)
            return 0;

        const umax log2x = cast(umax) find_last_set_bit(x);
        const umax log2b = cast(umax) find_last_set_bit(b);

        assert(log2x != 0);
        assert(log2b != 0);

        const umax lower_initial_estimate = log2x / (log2b + 1);
        const umax higher_initial_estimate = (log2x + log2b) / (log2b);

        umax curr_lo = lower_initial_estimate;
        umax curr_hi = higher_initial_estimate;

        while(true)
        {
            umax mid = (curr_lo + curr_hi) / 2;
            umax curr_approx = single_pow_by_squaring(b, mid);

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


    template <typename T>
    struct Slice
    {
        T* data = nullptr;
        size_t size = 0;

        func& operator[](size_t index) const noexcept { assert(index < size && "index out of range"); return data[index]; }
        func& operator[](size_t index) noexcept       { assert(index < size && "index out of range"); return data[index]; }

        constexpr operator Slice<const T>() const noexcept { return Slice<const T>{this->data, this->size};}
    };

    template <typename T>
    func slice_size(Slice<T> items, size_t from, size_t count)
    {
        assert(from <= items.size && from + count <= items.size && "sliced portion must be entirely within the given items");
        return Slice<T>{items.data + from, count};
    }

    template <typename T>
    func slice_range(Slice<T> items, size_t from, size_t to)
    {
        assert(to >= from);
        return slice_size<T>(items, from, to - from);
    }

    template <typename T>
    func slice(Slice<T> items, size_t from) -> Slice<T> 
    {
        return slice_size<T>(items, from, items.size - from);
    }

    template <typename T>
    func trim(Slice<T> items, size_t to_size) -> Slice<T> 
    {
        return slice_size<T>(items, 0, to_size);
    }

    template <typename T>
    func last(Slice<T> slice) -> T
    {
        return slice[slice.size - 1];
    }

    template <typename T>
    func last(Slice<T>* slice) -> T*
    {
        return &(*slice)[slice->size - 1];
    }
}
#endif // !TINY_NUM_FIRST_TIME_INCLUDED


#ifdef INCLUDED_NEW_TYPE

#ifndef TINY_NUM_LIB_NAME
#ifdef INCLUDED_TINY_NUM_TYPE_DEF
#define TINY_NUM_LIB_NAME tiny_num
#else
#error "name must be provided!"
#endif
#endif 

//so that the IDE displays something short...
#define tiny TINY_NUM_LIB_NAME
namespace tiny
{
    using namespace helpers;

    using Digit = INCLUDED_NEW_TYPE;

    using Slice = helpers::Slice<Digit>;
    using CSlice = helpers::Slice<const Digit>;

    struct Single_Overflow
    {
        Digit value = 0;
        Digit overflow = 0;

        bool constexpr operator ==(Single_Overflow const&) const noexcept = default;
    };

    struct Batch_Overflow
    {
        Slice slice;
        Digit overflow = 0;

        bool constexpr operator ==(Batch_Overflow const&) const noexcept = default;
    };

    func are_aliasing(CSlice left, CSlice right) -> bool
    { 
        uintptr_t left_pos = cast(uintptr_t) left.data;
        uintptr_t right_pos = cast(uintptr_t) right.data;
        if(right_pos < left_pos)
        {
            //[ right ]      [ left ]
            //[ ?? right size ?? ]

            uintptr_t diff = left_pos - right_pos;
            return diff < right.size;
        }
        else
        {
            //[ left ]      [ right ]
            //[ ?? left size ?? ]
            uintptr_t diff = right_pos - left_pos;
            return diff < left.size;
        }
    }

    func are_one_way_aliasing(CSlice before, CSlice after) -> bool
    { 
        return (before.data + before.size > after.data) && (after.data > before.data);
    }

    proc copy_slice(Slice* to, CSlice from) -> void
    {
        assert(to->size == from.size && "sizes must match");
        memmove(to->data, from.data, from.size * sizeof(Digit));
    }

    proc null_slice(Slice* to)
    {
        memset(to->data, 0, to->size * sizeof(Digit));
    }

    constexpr size_t DIGIT_BIT_SIZE = sizeof(Digit) * CHAR_BIT;
    constexpr size_t UMAX_BIT_SIZE = sizeof(umax) * CHAR_BIT;
    constexpr size_t DIGIT_HALF_BIT_SIZE = DIGIT_BIT_SIZE / 2;

    func set_bit(Digit field, size_t bit_pos, Digit val) -> Digit {
        assert(val == 0 || val == 1);
        assert(bit_pos < DIGIT_BIT_SIZE);

        const Digit bit = val << bit_pos;
        const Digit mask = cast(Digit) 1 << bit_pos;
        const Digit res = (field & ~mask) | bit;
        return res;
    };

    func get_bit(Digit field, size_t bit_pos) -> Digit 
    {
        return (field >> bit_pos) & 1u;
    };

    func set_nth_bit(Slice* num, size_t i, Digit val)
    {
        const size_t digit_i = i / DIGIT_BIT_SIZE;
        const size_t bit_i = i % DIGIT_BIT_SIZE;

        assert(digit_i < num->size);

        Digit* digit = &(*num)[digit_i];
        *digit = set_bit(*digit, bit_i, val);
    };

    func get_nth_bit(CSlice num, size_t i) -> Digit {
        const size_t digit_i = i / DIGIT_BIT_SIZE;
        const size_t bit_i = i % DIGIT_BIT_SIZE;

        return get_bit(num[digit_i], bit_i);
    };

    func high_mask(size_t index = DIGIT_HALF_BIT_SIZE) -> Digit {
        assert(index < DIGIT_BIT_SIZE);
        const Digit full_mask = -1;
        return cast(Digit) (full_mask << index);
    }

    func low_mask(size_t index = DIGIT_HALF_BIT_SIZE) -> Digit {
        return cast(Digit) ~high_mask(index);
    }

    func high_bits(Digit value, size_t index = DIGIT_HALF_BIT_SIZE) -> Digit {
        assert(index < DIGIT_BIT_SIZE);
        return cast(Digit) (value >> index);
    }

    func low_bits(Digit value, size_t index = DIGIT_HALF_BIT_SIZE) -> Digit {
        return cast(Digit) (value & low_mask(index));
    };  

    func combine_bits(Digit low, Digit high, size_t index = DIGIT_HALF_BIT_SIZE) -> Digit {
        assert(index < DIGIT_BIT_SIZE);
        return cast(Digit) (low_bits(low, index) | (high << index));
    };

    func dirty_combine_bits(Digit low, Digit high, size_t index = DIGIT_HALF_BIT_SIZE) -> Digit {
        assert(index < DIGIT_BIT_SIZE);
        assert(high_bits(low, index) == 0 && "low must not have high bits use combine_bits instead");
        return cast(Digit) (low | (high << index));
    };

    func find_last_set_digit(CSlice num) -> size_t
    {
        size_t i = num.size;
        for(; i-- > 0;)
            if(num[i] != 0)
                break;

        return i;
    }

    func find_first_set_digit(CSlice num) -> size_t
    {
        size_t i = 0;
        for(; i < num.size; i++)
            if(num[i] != 0)
                break;

        return i;
    }

    func striped_trailing_zeros(Slice num) -> Slice
    {
        return trim(num, find_last_set_digit(num) + 1);
    }

    func striped_trailing_zeros(CSlice num) -> CSlice
    {
        return trim(num, find_last_set_digit(num) + 1);
    }

    func is_striped_number(CSlice num) -> bool
    {
        if(num.size == 0)
            return true;

        return last(num) != 0;
    }

    func to_number(CSlice bignum) -> umax
    {
        CSlice stripped = striped_trailing_zeros(bignum);

        const size_t total_size = stripped.size * sizeof(Digit);
        assert(total_size <= sizeof(umax));

        umax out = 0;
        for(size_t i = 0; i < stripped.size; i++)
        {
            size_t shift_by = i * DIGIT_BIT_SIZE;
            assert(shift_by < UMAX_BIT_SIZE);
            Digit curr = stripped[i];
            out |= cast(umax) curr << shift_by;
        }

        return out;
    }

    proc from_number(Slice* bignum, umax from) -> Slice
    {
        constexpr size_t count = sizeof(umax) / sizeof(Digit);
        assert(bignum->size >= count);

        size_t i = 0;
        for(; i < count; i++)
        {
            umax curr_digit = from >> (i * DIGIT_BIT_SIZE);
            if(curr_digit == 0)
                break;

            (*bignum)[i] = cast(Digit) curr_digit;
        }

        return trim(*bignum, i);
    }

    func single_add_overflow_any(Digit carry, Digit a1, Digit a2, Digit a3) -> Single_Overflow
    {
        //Both add and substract work by splittin each 'digit' (type T) in half 
        // and performing desired operation then using the possible oveflow as a carry
        // for the next operation. 
        //
        //This efectively doubles the ammount of operations that need to be performed
        // but from my testing it is still faster than manually checking for overflow and
        // pathcing up (current approach requires no ifs in inner loops).
        //
        //We split each digit in exactly half since that guarantees that each half will be 
        // properly aligned to fit into a register 
        // (except T = u8 but those will not be used for more than testing purposes (I hope))

        assert(high_bits(carry) == 0);
        const Digit res_low = low_bits(a1) + low_bits(a2) + low_bits(a3) + carry;
        const Digit middle_carry = high_bits(res_low);

        const Digit res_high = high_bits(a1) + high_bits(a2) + high_bits(a3) + middle_carry;
        const Digit new_carry = high_bits(res_high);

        return Single_Overflow{combine_bits(res_low, res_high), new_carry};
    }

    func single_sub_overflow_any(Digit carry, Digit left, Digit a2, Digit a3) -> Single_Overflow
    {
        const Digit res_low = high_mask() + low_bits(left) - (low_bits(a2) + low_bits(a3)) - carry;
        const Digit middle_carry = low_mask() - high_bits(res_low);

        const Digit res_high = high_mask() + high_bits(left) - (high_bits(a2) + high_bits(a3)) - middle_carry;
        const Digit new_carry = low_mask() - high_bits(res_high);

        return Single_Overflow{combine_bits(res_low, res_high), new_carry};
    }

    func single_add_overflow(Digit left, Digit right, Digit carry = 0) -> Single_Overflow 
    {
        return single_add_overflow_any(carry, left, right, 0);
    }

    func single_sub_overflow(Digit left, Digit right, Digit carry = 0) -> Single_Overflow 
    {
        return single_sub_overflow_any(carry, left, right, 0);
    }

    func single_add_no_overflow(Digit a1, Digit a2 = 0, Digit a3 = 0, Digit a4 = 0, Digit a5 = 0) -> Digit
    {
        #ifdef NDEBUG
        return a1 + a2 + a3 + a4 + a5;
        #else
        const Digit arr[] = {a1, a2, a3, a4, a5};
        Digit sum = a1;

        for(size_t i = 1; i < 5; i++)
        {
            Digit new_sum = sum + arr[i];
            assert(new_sum >= sum && "must not overflow");
            sum = new_sum;
        }

        return sum;
        #endif // NDEBUG
    }

    func batch_add_or_sub_overflow_short(Slice* to, CSlice left, Digit carry, bool is_addition, Location location = Location::OUT_OF_PLACE, size_t from = 0) -> Batch_Overflow
    {
        assert(to->size >= left.size);
        assert(are_one_way_aliasing(left, *to) == false);

        if(location == Location::IN_PLACE)
            assert(to->data == left.data && "in place must be in place");

        const Slice trimmed_to = trim(*to, left.size);
        size_t j = from;
        for (; carry != 0 && j < left.size != 0; j++)
        {
            const Digit digit = left[j];
            bool carry_consumed = false;
            Digit patch_res = 0;
            if(is_addition)
            {
                patch_res = digit + carry;
                carry = !(patch_res > digit);
            }
            else
            {
                patch_res = digit - carry;
                carry = !(patch_res < digit);
            }

            trimmed_to[j] = patch_res;
        }

        if(location == Location::OUT_OF_PLACE)
        {
            Slice remainign_to = slice(trimmed_to, j);
            copy_slice(&remainign_to, slice(left, j));
        }

        return Batch_Overflow{trimmed_to, carry};
    }

    func batch_add_overflow_short(Slice* to, CSlice left, Digit right, Location location = Location::OUT_OF_PLACE, size_t from = 0) -> Batch_Overflow
    {
        return batch_add_or_sub_overflow_short(to, left, right, true, location, from);
    }

    func batch_sub_overflow_short(Slice* to, CSlice left, Digit right, Location location = Location::OUT_OF_PLACE, size_t from = 0) -> Batch_Overflow
    {
        return batch_add_or_sub_overflow_short(to, left, right, false, location, from);
    }

    func batch_add_overflow_long(Slice* to, CSlice left, CSlice right, Location location = Location::OUT_OF_PLACE, Digit carry_in = 0) -> Batch_Overflow
    {
        assert(to->size >= right.size);
        assert(to->size >= left.size);
        assert(are_one_way_aliasing(left, *to) == false);
        assert(high_bits(carry_in) == 0);

        if(left.size < right.size)
            swap(&left, &right);

        Digit carry = carry_in;
        for (size_t i = 0; i < right.size; i++)
        {
            Single_Overflow res = single_add_overflow(left[i], right[i], carry);
            (*to)[i] = res.value;
            carry = res.overflow;
        }

        return batch_add_overflow_short(to, left, carry, location, right.size);
    }

    func batch_sub_overflow_long(Slice* to, CSlice left, CSlice right, Location location = Location::OUT_OF_PLACE, Digit carry_in = 0) -> Batch_Overflow
    {
        assert(to->size >= right.size);
        assert(to->size >= left.size);
        assert(are_one_way_aliasing(left, *to) == false);
        assert(carry_in == 0 || carry_in == 1);

        Digit carry = carry_in;
        const size_t min_size = min(left.size, right.size);
        const size_t max_size = max(left.size, right.size);
        for (size_t i = 0; i < min_size; i++)
        {
            Single_Overflow oveflow = single_sub_overflow(left[i], right[i], carry);
            (*to)[i] = oveflow.value;
            carry = oveflow.overflow;
        }

        for (size_t i = left.size; i < right.size; i++)
        {
            Single_Overflow oveflow = single_sub_overflow(cast(Digit) 0, right[i], carry);
            (*to)[i] = oveflow.value;
            carry = oveflow.overflow;
        }

        if(right.size < left.size)
            return batch_sub_overflow_short(to, left, carry, location, right.size);

        return Batch_Overflow{trim(*to, max_size), carry};
    }

    func single_complement_overflow(Digit val, Digit carry) -> Single_Overflow
    {
        const Digit val_inv = ~val;

        const Digit complement_low = low_bits(val_inv) + carry;
        const Digit middle_carry = high_bits(complement_low);

        const Digit complement_high = high_bits(val_inv) + middle_carry;
        const Digit new_carry = high_bits(complement_high);

        return Single_Overflow{combine_bits(complement_low, complement_high), new_carry};
    }

    proc batch_complement_overflow(Slice* to, CSlice left, Digit carry_in = 1) -> Batch_Overflow
    {
        assert(to->size >= left.size);
        assert(is_striped_number(left));
        assert(are_one_way_aliasing(left, *to) == false);
        assert(carry_in == 1 || carry_in == 0);

        Digit carry = carry_in;
        for(size_t i = 0; i < left.size; i++)
        {
            Single_Overflow res = single_complement_overflow(left[i], carry);
            (*to)[i] = res.value;
            carry = res.overflow;
        }

        return Batch_Overflow{trim(*to, left.size), carry};
    }

    func single_shift_up_overflow(Digit low_item, size_t by_bits, Digit high_item) -> Single_Overflow
    {
        assert(by_bits < DIGIT_BIT_SIZE && by_bits != 0 && "shift must be in valid range (cannot be zero because its unimplmentable without if)");
        const size_t remaining_bits = DIGIT_BIT_SIZE - by_bits;

        // <--- loop iter direction <---
        // low_item | high_item
        // 0 1 0 1  | 1 0 0 1 >> by_bits: 3 => remaining_bits: 1
        // <-B-><R>   <-B-><R>
        //   <-L->   <H><-O->
        //
        // B:by_bits  R:remaining_bits
        // L:low      H:high
        // O:out (shifted_out)
        // 
        //            <-L-><H>   <-O->
        // 0 0 0 0 |  1 0 1 1 |  0 0 1 0
        //            <----->    <----->
        //             value       out

        const Digit low = high_bits(low_item, remaining_bits);
        const Digit high = low_bits(high_item, remaining_bits);

        const Digit composed = low | high << by_bits;

        const Digit out = high_bits(high_item, remaining_bits);
        const Digit shifted_out = out;
        return Single_Overflow{composed, shifted_out};
    }

    func single_shift_down_overflow(Digit low_item, size_t by_bits, Digit high_item) -> Single_Overflow
    {
        assert(by_bits < DIGIT_BIT_SIZE && by_bits != 0);
        const size_t remaining_bits = DIGIT_BIT_SIZE - by_bits;
        assert(remaining_bits < DIGIT_BIT_SIZE);

        // ---> loop iter direction --->
        // low_item | high_item
        // 0 1 0 1  | 1 0 0 1 << by_bits: 3 => remaining_bits: 1
        // <-B-><R>   <-B-><R>
        // <-O-><L>   <-H->
        //
        // B:by_bits  R:remaining_bits
        // L:low      H:high
        // O:out (shifted_out)
        // 
        //   <-O->   <L><-H->
        // 0 0 1 0 |  1 1 0 0  | 1 0 0 0
        // <----->    <----->
        //   out       value  

        
        const Digit low = high_bits(low_item, by_bits);
        const Digit high = low_bits(high_item, by_bits);

        const Digit composed = low | high << remaining_bits;

        const Digit out = low_bits(low_item, by_bits);
        const Digit shifted_out = out << remaining_bits;
        return Single_Overflow{composed, shifted_out};
    }

    //We allow shifting while iterating in both directions:
    // 
    // FOR SHIFT_UP: the FORWARD direction is useful in general case
    //  (for example when shifting array in place)
    // and the BACKWARD direction can be used while shifting by more then one element
    //  and writes to already passed memory are neccessary
    // 
    // FOR SHIFT_DOWN: its the polar opposite
    proc batch_shift_up_overflow(Slice* out, CSlice in, size_t by_bits, 
        Iter_Direction direction = Iter_Direction::FORWARD, Digit carry_in = 0) -> Batch_Overflow
    {
        assert(out->size >= in.size);
        assert(by_bits < DIGIT_BIT_SIZE);

        Slice out_trimmed = trim(*out, in.size);

        if(by_bits == 0)
        {
            copy_slice(&out_trimmed, in);
            return Batch_Overflow{out_trimmed, carry_in};
        }

        if(in.size == 0)
            return Batch_Overflow{out_trimmed, carry_in};

        if(direction == Iter_Direction::FORWARD)
        {
            assert(are_one_way_aliasing(in, *out) == false);
            Digit prev = carry_in;
            for (size_t i = 0; i < in.size; i++)
            {
                Digit curr = in[i];
                Single_Overflow res = single_shift_up_overflow(prev, by_bits, curr);
                (*out)[i] = res.value;
                prev = curr;
            }

            Single_Overflow carry_out = single_shift_up_overflow(prev, by_bits, 0);
            return Batch_Overflow{out_trimmed, carry_out.value};
        }
        else
        {
            assert(are_one_way_aliasing(*out, in) == false);
            const Digit shifted_out = in[in.size - 1];
            Digit prev = shifted_out;
            for (size_t i = in.size; i-- > 1; )
            {
                Digit curr = in[i - 1];
                (*out)[i] = single_shift_up_overflow(curr, by_bits, prev).value;
                prev = curr;
            }

            {
                Single_Overflow res = single_shift_up_overflow(carry_in, by_bits, prev);
                (*out)[0] = res.value;
            }


            Single_Overflow carry_out = single_shift_up_overflow(shifted_out, by_bits, 0);
            return Batch_Overflow{out_trimmed, carry_out.value};
        }
    }

    proc batch_shift_down_overflow(Slice* out, CSlice in, size_t by_bits, 
        Iter_Direction direction = Iter_Direction::FORWARD, Digit carry_in = 0) -> Batch_Overflow
    {
        assert(out->size >= in.size);
        assert(by_bits < DIGIT_BIT_SIZE);

        Slice out_trimmed = trim(*out, in.size);
        if(by_bits == 0)
        {
            //@TODO: Handle carry_in! (consistently with the other edge case)
            copy_slice(&out_trimmed, in);
            return Batch_Overflow{out_trimmed, carry_in};
        }

        if(in.size == 0)
            return Batch_Overflow{out_trimmed, carry_in};

        if(direction == Iter_Direction::FORWARD)
        {
            assert(are_one_way_aliasing(in, *out) == false);

            const Digit shifted_out = in[0];
            Digit prev = shifted_out;
            for (size_t i = 0; i < in.size - 1; i++)
            {
                Digit curr = in[i + 1];
                Single_Overflow res = single_shift_down_overflow(prev, by_bits, curr);
                (*out)[i] = res.value;
                prev = curr;
            }

            {
                Single_Overflow res = single_shift_down_overflow(prev, by_bits, carry_in);
                (*out)[in.size - 1] = res.value;
            }

            Single_Overflow carry_out = single_shift_down_overflow(0, by_bits, shifted_out);
            return Batch_Overflow{out_trimmed, carry_out.value};
        }
        else
        {
            assert(are_one_way_aliasing(*out, in) == false);

            Digit prev = carry_in;
            for (size_t i = in.size; i-- > 0;)
            {
                Digit curr = in[i];
                (*out)[i] = single_shift_down_overflow(curr, by_bits, prev).value;
                prev = curr;
            }

            Single_Overflow carry_out = single_shift_down_overflow(0, by_bits, prev);
            return Batch_Overflow{out_trimmed, carry_out.value};
        }
    }

    //return -1 if left < right
    //        0 if left == right
    //        1 if left > right   
    func compare(CSlice left, CSlice right) -> int
    {
        assert(is_striped_number(left));
        assert(is_striped_number(right));
        if(left.size < right.size)
            return -1;

        if(left.size > right.size)
            return 1;

        for(size_t i = left.size; i-- > 0;)
        {
            if(left[i] < right[i])
                return -1;
            if(left[i] > right[i])
                return 1;
        }

        return 0;
    }

    func is_equal(CSlice left, CSlice right) -> bool
    {
        if(left.size != right.size)
            return false;

        for(size_t i = 0; i < left.size; i++)
            if(left[i] != right[i])
                return false;

        return true;
    }

    enum class Mul_Overflow_Optims
    {
        NONE,
        HIGH_BITS_ONLY,
        LOW_BITS_ONLY,
    };

    func single_mul_overflow(Digit left, Digit right, Digit last_value, Mul_Overflow_Optims const& mul_optims = Mul_Overflow_Optims::NONE) -> Single_Overflow
    {
        //we do the oveflow multiplication by multiplying each digit normally and summing the overflow 
        // => 25  
        //    *5
        //  ----
        //    25 <- here 2 is 'carry over' to the next iteration 
        //   10  
        //  ----
        //   125 <- normally we add all partial results all at once but during our algorithm we do it single pass
        // 
        // Since we cant have only finite machine integers we think of each value as having two digits (slots)
        //  labeled low and high (l1 and l2 respectively) 
        // 
        // l = l1 + l2 * b   (left)
        // r = r1 + r2 * b   (right)
        // h... half bits
        // b... base ie 2^h => x * b == x << h
        //
        // l * r = (l1 + l2*b) * (r1 + r2*b)
        //       = l1r1 + l1r2*b + l2r1*b + l2r2*b^2 
        // 
        // Slots:  [1..2]  [2..3]   [2..3]  [3 .. 4]
        //         [1..2]  [2    ...    3]  [3 .. 4]
        //          s_cur     s_mixed        s_next
        // 
        // => s_mixed slot is half in the next iteration => we dont shift (cause we would shift out of type)
        // => s_next slot is entirely in the next iteration => also no shift
        // => even though each of the two individual s_mixed summants will not overflow (half multiply never overflows) the 
        //     sum might - then it is handled as next digit and added to s_next
        // 
        // s_mixed = s_mixed_ns * b   //ns means 'no shift'
        // s_next = s_next_ns * b 
        // 
        // => last_value can also overflow the current value and it too needs to be added to next
        
        size_t h = DIGIT_HALF_BIT_SIZE;

        Digit l = left;
        Digit r = right;

        Digit l1 = low_bits(l); 
        Digit l2 = high_bits(l); 
        Digit r1 = low_bits(r);  
        Digit r2 = high_bits(r); 

        //calculate constants
        Digit s_cur = 0;
        Digit s_mixed_ns = 0;
        Digit s_mixed_ns_overflow = 0;
        Digit s_next_ns = 0;

        //if we have either only low bits or high bits we dont have to do
        // one carry addition and one multiplication. This is pretty big difference when dealing with really large numbers
        // (there are also additional parts of the code that can be optimized away - such as shifting of s_mixed_ns_overflow etc...)
        if(mul_optims == Mul_Overflow_Optims::LOW_BITS_ONLY)
        {
            assert(high_bits(right) == 0);
            s_cur = l1*r1;
            s_mixed_ns = l2*r1;
        }
        else if(mul_optims == Mul_Overflow_Optims::HIGH_BITS_ONLY)
        {
            assert(low_bits(right) == 0);
            s_mixed_ns = l1*r2;
            s_next_ns = l2*r2;
        }
        else
        {
            Single_Overflow s_mixed_ns_res = single_add_overflow(l1*r2, l2*r1);

            s_mixed_ns = s_mixed_ns_res.value;
            s_mixed_ns_overflow = s_mixed_ns_res.overflow;
            s_cur = l1*r1;
            s_next_ns = l2*r2;
        }

        //split mixed
        assert(h < DIGIT_BIT_SIZE);

        Digit low_mixed = low_bits(s_mixed_ns) << h; // even though we are taking lower half of s_mixed we have to insert it to upper half
        Digit high_mixed = high_bits(s_mixed_ns) + (s_mixed_ns_overflow << h); //add overflow as another digit

        //distribute mixed to appropriate halves

        Single_Overflow curr_value = single_add_overflow_any(0, s_cur, low_mixed, last_value); //also add last_value to save ops
        Digit next_value = single_add_no_overflow(s_next_ns, high_mixed, curr_value.overflow); //also add the overflow

        return Single_Overflow{curr_value.value, next_value};
    }

    func single_is_power_of_two(Digit num)
    {
        return num != 0 && (num & (num - 1)) == 0;
    }

    proc batch_mul_overflow(Slice* to, CSlice left, Digit right, Optims const& optims, Digit carry_in = 0) -> Batch_Overflow
    {
        assert(to->size >= left.size);
        assert(are_one_way_aliasing(left, *to) == false);

        Digit carry = carry_in;
        Slice trimmed_to = trim(*to, left.size);

        if(optims.mul_consts)
        {
            if(left.size == 0)
                return Batch_Overflow{trimmed_to, carry_in};

            if(right == 0) 
                return Batch_Overflow{trim(*to, 0), 0};

            if(right == 1)
            {
                copy_slice(&trimmed_to, left);
                return Batch_Overflow{trimmed_to, 0};
            }
        }

        if(optims.mul_shift)
        {
            if(single_is_power_of_two(right))
            {
                //@NOTE: here carry is handled differently the in the rest of the mul algorhitm
                //  this doesnt affect how the algorhitm appears from the outside not does it affect the
                //  ability for batches to be chained together with carries but it forbids any kind of dependence
                //  on the carry in/out in this case. This is not ideal and should be changed so that the
                //  shift itself handles carry in similar way to mul
                size_t shift_bits = find_last_set_bit(right);
                return batch_shift_up_overflow(to, left, shift_bits, Iter_Direction::FORWARD, carry_in);
            }
        }

        #define do_carry_loop(operation)            \
        {                                           \
            for (size_t i = 0; i < left.size; i++)  \
            {                                       \
                const Single_Overflow res = operation;  \
                (*to)[i] = res.value;               \
                carry = res.overflow;               \
            }                                       \
        }                                           \

        if(optims.mul_half_bits && high_bits(right) == 0)
            do_carry_loop((single_mul_overflow(left[i], right, carry, Mul_Overflow_Optims::LOW_BITS_ONLY)))
        else if(optims.mul_half_bits && low_bits(right) == 0)
            do_carry_loop((single_mul_overflow(left[i], right, carry, Mul_Overflow_Optims::HIGH_BITS_ONLY)))
        else
            do_carry_loop((single_mul_overflow(left[i], right, carry)))

        #undef do_carry_loop
        return Batch_Overflow{trimmed_to, carry};
    }

    func single_div_overflow(Digit left, Digit right, Digit carry_in = 0) -> Single_Overflow
    {
        //The algorhitm works as follows (only in different base - we use base 10 for demosntartion)
        // 61 / 5 == 10 + 11 / 5 == 10 + 2 == 12
        // 
        //To go through this instinctive algorhitm we did the followig:
        // {161/5} == [16/5]*10 + (16%5)*10 + {1/5} = 3*10 + {(10 + 1) / 5} == 30 + {11/5} == ... == 32 
        // where {} means normal (ie infinitely precise fraction) used for yet uncalculated division

        //this however only works for carry and right thats smaller than half bits 
        // (else it doesnt get carried through properly through the modulos
        //  so for example single_div_overflow(0, 0xFF, 0xF0) should return 0xF0 but returns less) => assert
        assert(high_bits(carry_in) == 0 && "carry must be single digit");
        assert(high_bits(right) == 0 && "right must be single digit");
        assert(right != 0 && "cannot divide by zero");

        const Digit operand_high = combine_bits(high_bits(left), carry_in);
        const Digit res_high = operand_high / right;
        const Digit middle_carry = operand_high % right;

        const Digit operand_low = combine_bits(low_bits(left), middle_carry);
        const Digit res_low = operand_low / right;
        const Digit out_carry = operand_low % right;

        const Digit res = dirty_combine_bits(res_low, res_high);
        return Single_Overflow{res, out_carry};
    }

    proc batch_div_overflow(Slice* to, CSlice left, Digit right, Optims const& optims, Digit carry_in = 0) -> Batch_Overflow
    {
        assert(to->size >= left.size);
        assert(high_bits(right) == 0 && "only works for divisors under half bit size");
        assert(right != 0 && "cannot divide by zero");
        assert(are_one_way_aliasing(*to, left) == false);

        Slice trimmed_to = trim(*to, left.size);
        if(optims.div_consts)
        {
            if(right == 1)
            {
                copy_slice(&trimmed_to, left);
                return Batch_Overflow{trimmed_to, 0}; 
            }
        }

        if(optims.div_shift)
        {
            if(single_is_power_of_two(right))
            {
                size_t shift_bits = find_last_set_bit(right);
                Batch_Overflow shift_res = batch_shift_down_overflow(to, left, shift_bits, Iter_Direction::BACKWARD, carry_in);
                //because of the way shifting down result is defined we have to process the reuslt
                const Digit carry = shift_res.overflow >> (DIGIT_BIT_SIZE - shift_bits);
                return Batch_Overflow{trimmed_to, carry};
            }
        }

        Digit carry = carry_in;
        for (size_t i = left.size; i-- > 0; )
        {
            Single_Overflow res = single_div_overflow(left[i], right, carry);
            (*to)[i] = res.value;
            carry = res.overflow;
        }

        return Batch_Overflow{trimmed_to, carry};
    }

    proc batch_rem_overflow(CSlice left, Digit right, Optims const& optims, Digit carry_in = 0) -> Digit
    {
        assert(high_bits(right) == 0 && "only works for divisors under half bit size");
        assert(right != 0 && "cannot divide by zero");

        if(optims.rem_optims)
        {
            if(single_is_power_of_two(right))
            {
                size_t shift_bits = find_last_set_bit(right);
                if(left.size == 0)
                    return carry_in;

                return low_bits(left[0], shift_bits);
            }
        }

        Digit carry = carry_in;
        for (size_t i = left.size; i-- > 0; )
        {
            Single_Overflow res = single_div_overflow(left[i], right, carry);
            carry = res.overflow;
        }

        return carry;
    }

    func required_add_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return max(left_size, right_size) + 1;
    }

    proc add(Slice* to, CSlice left, CSlice right, Location location = Location::OUT_OF_PLACE)
    {
        assert(is_striped_number(left));
        assert(is_striped_number(right));
        assert(to->size >= required_add_out_size(left.size, right.size));

        Batch_Overflow res = batch_add_overflow_long(to, left, right, location);
        if(res.overflow == 0)
        {
            if(location == Location::IN_PLACE)
                return striped_trailing_zeros(res.slice);
            
            assert(is_striped_number(res.slice));
            return res.slice;
        }

        size_t slice_size = res.slice.size;
        Slice trimmed_to = trim(*to, slice_size + 1);
        trimmed_to[slice_size] = res.overflow;

        assert(is_striped_number(trimmed_to));

        return trimmed_to;
    }

    func required_sub_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return left_size;
    }

    proc sub(Slice* to, CSlice left, CSlice right, Location location = Location::OUT_OF_PLACE)
    {
        assert(is_striped_number(left));
        assert(is_striped_number(right));
        assert(to->size >= required_sub_out_size(left.size, right.size));

        Batch_Overflow res = batch_sub_overflow_long(to, left, right, location);
        assert(res.overflow == 0 && "left must be bigger than right");

        return striped_trailing_zeros(res.slice);
    }

    func required_mul_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return left_size + right_size;
    }

    proc mul(Slice* to, CSlice left, Digit right, Optims const& optims)
    {
        assert(is_striped_number(left));
        assert(to->size >= required_sub_out_size(left.size, 1));

        Batch_Overflow res = batch_mul_overflow(to, left, right, optims);
        if(res.overflow == 0)
        {
            assert(is_striped_number(res.slice));
            return res.slice;
        }

        size_t slice_size = res.slice.size;
        Slice trimmed_to = trim(*to, slice_size + 1);
        trimmed_to[slice_size] = res.overflow;

        assert(is_striped_number(trimmed_to));
        return trimmed_to;
    }

    func required_div_quotient_size(size_t num_size, size_t den_size) -> size_t
    {
        if(num_size < den_size)
            return 0;

        return num_size - den_size + 1;
    }

    func required_div_remainder_size(size_t num_size, size_t den_size) -> size_t
    {
        return min(den_size + 1, num_size);;
    }

    struct Div_Result
    {
        Slice quotient;
        Slice remainder;
    };

    proc div_bit_by_bit(Slice* quotient, Slice* remainder, CSlice num, CSlice den, Optims const& optims, bool DO_QUOTIENT = true) -> Div_Result
    {
        assert(is_striped_number(num));
        assert(is_striped_number(den));
        assert(are_aliasing(*quotient, num) == false);
        assert(are_aliasing(*remainder, num) == false);
        assert(are_aliasing(*quotient, den) == false);
        assert(are_aliasing(*remainder, den) == false);

        const size_t required_quot_size = DO_QUOTIENT ? required_div_quotient_size(num.size, den.size) : 0;
        const size_t required_rem_size = required_div_remainder_size(num.size, den.size);

        if(DO_QUOTIENT)
            assert(are_aliasing(*quotient, *remainder) == false);

        assert(den.size > 0 && "cannot divide by zero!");

        if(num.size == 0)
        {
            Slice stripped_quotient = trim(*quotient, 0);
            Slice stripped_remainder = trim(*remainder, 0);
            return Div_Result{stripped_quotient, stripped_remainder};
        }

        if(num.size < den.size)
        {
            assert(remainder->size >= num.size);
            Slice stripped_remainder = trim(*remainder, num.size);
            Slice stripped_quotient = trim(*quotient, 0);

            copy_slice(&stripped_remainder, num);
            return Div_Result{stripped_quotient, stripped_remainder};
        }

        Slice trimmed_remainder = trim(*remainder, required_rem_size);
        Slice trimmed_quotient = trim(*quotient, required_quot_size);    

        if(DO_QUOTIENT)
            null_slice(&trimmed_quotient);

        null_slice(&trimmed_remainder);
        if(den.size == 1 && high_bits(last(den)) == 0)
        {
            assert(trimmed_remainder.size > 0 && "at this point should not be 0");

            Slice stripped_quotient;
            Digit remainder_val;
            if(DO_QUOTIENT)
            {
                Batch_Overflow div_res = batch_div_overflow(&trimmed_quotient, num, last(den), optims);
                stripped_quotient = striped_trailing_zeros(div_res.slice);
                remainder_val = div_res.overflow;
            }
            else
            {
                remainder_val = batch_rem_overflow(num, last(den), optims);
                stripped_quotient = trim(trimmed_quotient, 0);
            }

            size_t remainder_size = 0;
            if(remainder_val != 0)
            {
                trimmed_remainder[0] = remainder_val;
                remainder_size = 1;
            }

            Slice stripped_remainder = trim(trimmed_remainder, remainder_size);
            return Div_Result{stripped_quotient, stripped_remainder};
        }

        Slice curr_remainder = trim(trimmed_remainder, 0);
        for(size_t i = DIGIT_BIT_SIZE * num.size; i-- > 0; ) 
        {
            size_t last_size = curr_remainder.size;
            Batch_Overflow shift_res = batch_shift_up_overflow(&curr_remainder, curr_remainder, 1);
            const Digit num_ith_bit = get_nth_bit(num, i);

            if(shift_res.overflow != 0 || (last_size == 0 && num_ith_bit == 1))
            {
                curr_remainder = trim(trimmed_remainder, last_size + 1);
                curr_remainder[last_size] = shift_res.overflow;
            }

            set_nth_bit(&trimmed_remainder, 0, num_ith_bit);
            if(compare(curr_remainder, den) >= 0)
            {
                Batch_Overflow sub_res = batch_sub_overflow_long(&curr_remainder, curr_remainder, den, Location::IN_PLACE);
                curr_remainder = striped_trailing_zeros(curr_remainder);
                assert(sub_res.overflow == 0 && "should not overflow");
                    
                if(DO_QUOTIENT)
                    set_nth_bit(&trimmed_quotient, i, 1);
            }
        }

        assert(is_striped_number(curr_remainder));
        Slice stripped_quotient = DO_QUOTIENT 
            ? striped_trailing_zeros(trimmed_quotient)
            : trim(trimmed_quotient, 0);

        return Div_Result{stripped_quotient, curr_remainder};
    }

    proc div(Slice* quotient, Slice* remainder, CSlice num, CSlice den, Optims const& optims) -> Div_Result
    {
        return div_bit_by_bit(quotient, remainder, num, den, optims);
    }

    proc rem(Slice* remainder, CSlice num, CSlice den, Optims const& optims) -> Slice
    {
        Slice quotient = {nullptr, 0};
        Div_Result div_res = div_bit_by_bit(&quotient, remainder, num, den, optims, false);

        assert(div_res.quotient.size == 0 && "in case of only remainder the quotient size should be 0");
        assert(div_res.quotient.data == nullptr && "and data nullptr as we set it");
        return div_res.remainder;
    }

    proc batch_fused_mul_add_overflow(Slice* to, CSlice added, CSlice multiplied, Digit coeficient, 
        Optims const& optims, const Location location = Location::OUT_OF_PLACE, Digit add_carry = 0, Digit mul_carry = 0) -> Batch_Overflow
    {   
        assert(added.size >= multiplied.size);
        assert(to->size >= added.size);
        assert(are_one_way_aliasing(added, *to) == false);
        assert(are_one_way_aliasing(multiplied, *to) == false);

        if(optims.mul_consts)
        {
            Slice trimmed_to = trim(*to, added.size);
            if(coeficient == 0) 
            {
                if (location == Location::IN_PLACE)
                    return Batch_Overflow{trim(*to, 0), 0};
                else
                {
                    copy_slice(&trimmed_to, added);
                    return Batch_Overflow{trimmed_to, 0};
                }
            }

            if(coeficient == 1)
                return batch_add_overflow_long(to, added, multiplied, location);
        }

        Slice trimmed_to = trim(*to, added.size);
        if(optims.mul_shift 
            && optims.mul_consts //zero check must be performed (else shifting is UB!)
            && single_is_power_of_two(coeficient))
        {
            size_t shift_bits = find_last_set_bit(coeficient);
            const Slice trimmed_to = trim(*to, added.size);
            Digit prev_muled = mul_carry;
            for (size_t j = 0; j < multiplied.size; j++)
            {
                const Digit curr_muled = multiplied[j];
                const Digit curr_added = added[j];

                Single_Overflow mul_res = single_shift_up_overflow(prev_muled, shift_bits, curr_muled);
                Single_Overflow add_res = single_add_overflow(curr_added, mul_res.value, add_carry);

                trimmed_to[j] = add_res.value;

                prev_muled = curr_muled;
                add_carry = add_res.overflow;
            }

            mul_carry = single_shift_up_overflow(prev_muled, shift_bits, 0).value;
        }
        else
        {
            for (size_t j = 0; j < multiplied.size; j++)
            {
                Digit curr_muled = multiplied[j];
                Digit curr_added = added[j];

                Single_Overflow mul_res = single_mul_overflow(curr_muled, coeficient, mul_carry);
                Single_Overflow add_res = single_add_overflow(curr_added, mul_res.value, add_carry);

                trimmed_to[j] = add_res.value;

                mul_carry = mul_res.overflow;
                add_carry = add_res.overflow;
            }
        }
        
        const Digit combined_carry = single_add_no_overflow(mul_carry, add_carry);
        return batch_add_overflow_short(&trimmed_to, added, combined_carry, location, multiplied.size);
    }

    proc mul(Slice* to, Slice* aux, CSlice left, CSlice right, Optims const& optims, size_t depth = 0) -> Slice;

    func required_mul_quadratic_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return required_mul_out_size(left_size, right_size);
    }

    func required_mul_quadratic_aux_size(size_t left_size, size_t right_size) -> size_t
    {
        return required_mul_out_size(left_size, right_size);
    }

    proc mul_quadratic(Slice* to, Slice* temp, CSlice left, CSlice right, Optims const& optims) -> Slice
    {
        assert(are_aliasing(*to, left) == false);
        assert(are_aliasing(*to, right) == false);
        assert(are_aliasing(*temp, left) == false);
        assert(are_aliasing(*temp, right) == false);
        assert(are_aliasing(*to, *temp) == false);

        if(left.size < right.size)
            swap(&left, &right);

        const size_t max_result_size = required_mul_out_size(left.size, right.size);
        const size_t max_temp_size = required_mul_out_size(left.size, right.size);
        assert(to->size >= max_result_size);
        assert(temp->size >= max_temp_size);
        
        Slice trimmed_to = trim(*to, max_result_size);
        Slice trimmed_temp = trim(*temp, max_temp_size);
        null_slice(&trimmed_to);

        for(size_t i = 0; i < right.size; i++)
        {
            Batch_Overflow mul_res = batch_mul_overflow(&trimmed_temp, left, right[i], optims);

            //batch_mul_overflow is limited to left.size elems => in case of overflow
            // even if the overflow would fit into temp it is not added => we add it back
            size_t mul_size = mul_res.slice.size;
            Slice mul_slice = trim(trimmed_temp, mul_size + 1); 
            mul_slice[mul_size] = mul_res.overflow;

            Slice shifted_to = slice(trimmed_to, i);
            assert(shifted_to.size >= mul_slice.size);

            Batch_Overflow add_res = batch_add_overflow_long(&shifted_to, shifted_to, mul_slice, Location::IN_PLACE);
            assert(add_res.overflow == 0 && "no final carry should be left");
        }

        return striped_trailing_zeros(trimmed_to);
    }

    func required_mul_quadratic_fused_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return required_mul_out_size(left_size, right_size);
    }

    proc mul_quadratic_fused(Slice* to, CSlice left, CSlice right, Optims const& optims) -> Slice
    {
        const size_t max_result_size = required_mul_out_size(left.size, right.size);
        assert(are_aliasing(*to, left) == false);
        assert(are_aliasing(*to, right) == false);
        assert(to->size >= max_result_size);

        if(left.size < right.size)
            swap(&left, &right);


        Slice trimmed_to = trim(*to, max_result_size);

        //we null only the place for the first multiplication / addition
        // the rest is nulled implicitly through the overflow handling code
        // this saves us one additinal iteration through the data
        Slice first_iter = trim(*to, left.size);
        null_slice(&first_iter);

        for(size_t i = 0; i < right.size; i++)
        {
            Slice shifted_to = slice_size(trimmed_to, i, left.size);

            Batch_Overflow mul_add_res = batch_fused_mul_add_overflow(&shifted_to, shifted_to, left, right[i], optims, Location::IN_PLACE);
            trimmed_to[i + left.size] = mul_add_res.overflow;
        }

        return striped_trailing_zeros(trimmed_to);
    }

    func required_mul_karatsuba_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return required_mul_out_size(left_size, right_size);
    }

    func required_mul_karatsuba_aux_size(size_t left_size, size_t right_size, size_t recursion_depth) -> size_t
    {
        if(left_size < right_size)
            swap(&left_size, &right_size);

        size_t base = max(left_size, right_size) / 2; 
        size_t capped_digits = min(right_size, base);

        size_t x1 = base;
        size_t x2 = left_size - base;

        size_t y1 = capped_digits;
        size_t y2 = right_size - capped_digits;

        size_t required_add_x_sum_size = required_add_out_size(x1, x2);
        size_t required_add_y_sum_size = required_add_out_size(y1, y2);

        size_t required_z2_size = -1;
        if(recursion_depth == 0)
            required_z2_size = required_mul_out_size(required_add_x_sum_size, required_add_y_sum_size);
        else
            required_z2_size = required_mul_karatsuba_aux_size(required_add_x_sum_size, required_add_y_sum_size, recursion_depth - 1);

        return required_z2_size 
            + required_add_x_sum_size 
            + required_add_y_sum_size;
    }

    proc mul_karatsuba(Slice* to, Slice* aux, CSlice left, CSlice right, Optims const& optims, size_t depth = 0, bool is_run_alone = true) -> Slice
    {
        assert(is_striped_number(left));
        assert(is_striped_number(right));

        assert(are_aliasing(*to, left) == false);
        assert(are_aliasing(*aux, left) == false);
        assert(are_aliasing(*to, right) == false);
        assert(are_aliasing(*aux, right) == false);
        assert(are_aliasing(*aux, *to) == false);

        CSlice x = left;
        CSlice y = right;

        //Perfomrs the karatsuba multiplicaton step
        //x*y = (x1*b + x2)(y1*b + y2) = (x1*y1)b^2 + (x1*y2 + x2*y1)b + (x2*y2)
        //    = z1*b^2 + z2*b + z3

        // z1 = x1*y1
        // z3 = x2*y2
        // z2 = x1*y2 + x2*y1 = (x1 + x2)(y1 + y2) - z1 - z3

        if(x.size < y.size)
            swap(&x, &y);

        if(is_run_alone)
        {
            if(x.size == 0 || y.size == 0)
                return trim(*to, 0);

            if(y.size == 1)
                return mul(to, x, y[0], optims);
        }

        //size_t split_digits = min(x.size / 2, y.size / 2);
        //we split according to the bigger one. This might seem odd but in the case where we will split too much
        // and y1 will be 0 we will have one less multiplication to worry about
        size_t base = max(x.size, y.size) / 2; 
        size_t capped_digits = min(y.size, base);
        assert(base != 0);

        CSlice x1 = slice(x, base);
        CSlice x2 = striped_trailing_zeros(trim(x, base));

        CSlice y1 = slice(y, capped_digits);
        CSlice y2 = striped_trailing_zeros(trim(y, capped_digits));


        size_t required_out_size = required_mul_out_size(x.size, y.size);
        Slice trimmed_out = trim(*to, required_out_size);

        size_t required_add_x_sum_size = required_add_out_size(x1.size, x2.size);
        size_t required_add_y_sum_size = required_add_out_size(y1.size, y2.size);

        size_t required_z1_size = required_mul_out_size(x1.size, y1.size);
        size_t required_z2_size = required_mul_out_size(required_add_x_sum_size, required_add_y_sum_size); //after multiplies will be even less
        size_t required_z3_size = required_mul_out_size(x2.size, y2.size);

        size_t total_used_aux_size = required_add_x_sum_size + required_add_y_sum_size + required_z2_size;
        size_t returned_required_aux = required_mul_karatsuba_aux_size(left.size, right.size, 0);
        assert(total_used_aux_size <= returned_required_aux);

        Slice remaining_aux = slice(*aux, total_used_aux_size);

        size_t z1_from_index = trimmed_out.size - required_z1_size;

        //we place z1 and z3 directly into the output buffer 
        // into appropiate positions as if multiplied by base
        Slice z1_slot = slice(trimmed_out, z1_from_index);
        Slice z3_slot = trim(trimmed_out, required_z3_size);

        assert(are_aliasing(z1_slot, z3_slot) == false);

        //we multyply into z3 and null the area between it and begginign of z1 
        // (ie the place where z2 will get added to)
        Slice z3 = mul(&z3_slot, &remaining_aux, x2, y2, optims, depth + 1);
        Slice z3_to_z1 = slice_range(trimmed_out, z3.size, z1_from_index);
        null_slice(&z3_to_z1);

        //we multiply into z1 and null the remaining size to the end of trimmed_out
        // this way trimmed out should look like the following:
        // 000[ z1 ]000[ z3 ]
        Slice z1 = mul(&z1_slot, &remaining_aux, x1, y1, optims, depth + 1);
        Slice z1_up = slice(trimmed_out, z1.size + z1_from_index);
        null_slice(&z1_up);

        size_t used_to = 0;
        Slice z2_slot = slice_size(*aux, used_to, required_z2_size);
        used_to += required_z2_size;

        Slice x_sum_slot = slice_size(*aux, used_to, required_add_x_sum_size);
        used_to += required_add_x_sum_size;

        Slice y_sum_slot = slice_size(*aux, used_to, required_add_y_sum_size);
        used_to += required_add_y_sum_size;

        Slice x_sum = add(&x_sum_slot, x1, x2);
        Slice y_sum = add(&y_sum_slot, y1, y2);

        Slice x_sum_y_sum = mul(&z2_slot, &remaining_aux, x_sum, y_sum, optims, depth + 1);
        Slice x_sum_y_sum_m_z1 = sub(&x_sum_y_sum, x_sum_y_sum, z1, Location::IN_PLACE);
        Slice z2 = sub(&x_sum_y_sum_m_z1, x_sum_y_sum_m_z1, z3, Location::IN_PLACE);

        //instead of multiplying z2 by base we add it to the appropriate position
        Slice out_z2_up = slice(trimmed_out, base);

        Batch_Overflow add_res = batch_add_overflow_long(&out_z2_up, out_z2_up, z2, Location::IN_PLACE);
        assert(add_res.overflow == 0 && "should not overflow");

        Slice out = striped_trailing_zeros(trimmed_out);

        return out;
    }

    proc mul(Slice* to, Slice* aux, CSlice left, CSlice right, Optims const& optims, size_t depth) -> Slice
    {
        assert(is_striped_number(left));
        assert(is_striped_number(right));
        assert(are_aliasing(*to, left) == false);
        assert(are_aliasing(*to, right) == false);
        assert(to->size >= required_mul_out_size(left.size, right.size));

        Slice ret;
    
        size_t max_size = -1;
        size_t min_size = -1;

        if(left.size > right.size)
        {
            max_size = left.size;
            min_size = right.size;
        }
        else
        {
            max_size = left.size;
            min_size = right.size;
        }

        if(min_size >= optims.mul_quadratic_single_below_size 
            && max_size >= optims.mul_quadratic_both_below_size 
            && depth < optims.max_recursion_depth)
        {
            size_t required_aux = required_mul_karatsuba_aux_size(max_size, min_size, 1);
            if(aux->size >= required_aux)
            {
                bool is_run_alone = optims.mul_quadratic_single_below_size <= 1;
                ret = mul_karatsuba(to, aux, left, right, optims, depth, is_run_alone);
                assert(is_striped_number(ret));
                return ret;
            }
        }
        
        ret = mul_quadratic_fused(to, left, right, optims);
        assert(is_striped_number(ret));
        return ret;
    }

    func log2(CSlice left) -> size_t
    {
        assert(is_striped_number(left));
        if(left.size == 0)
            return 0;

        size_t last_log = find_last_set_bit(last(left));
        return last_log + (left.size - 1) * DIGIT_BIT_SIZE;
    }

    func required_pow_out_size(size_t num_bit_size, umax power, size_t single_digit_bit_size) -> size_t
    {
        if(power == 0 || num_bit_size == 0)
            return 1;

        size_t bit_size = (num_bit_size + 1) * power;
        size_t item_size = div_round_up(bit_size, single_digit_bit_size) + 1;

        return item_size;
    }

    func required_pow_by_squaring_single_aux_swap_size(size_t num_bit_size, umax power, size_t single_digit_bit_size) -> size_t
    {
        if(power == 0)
            return 0;

        size_t iters = find_last_set_bit(power);
        size_t bits = (num_bit_size + 1);
        //we square the number on every iteration => the number of bits doubles => 2^iters == 1 << iters
        size_t bit_size = bits << iters; // == bits * 2^iters == bits * 2^[log2(power)] ~~ bits*power
        size_t item_size = div_round_up(bit_size, single_digit_bit_size) + 1;
        return item_size;
    }

    func required_pow_by_squaring_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size) -> size_t
    {
        if(power == 0)
            return 0;
        
        size_t square = required_pow_by_squaring_single_aux_swap_size(num_bit_size, power, single_digit_bit_size);
        size_t out_swap = required_pow_out_size(num_bit_size, power, single_digit_bit_size);

        return square * 2 + out_swap;
    }

    func optimal_pow_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size) {
        return required_pow_by_squaring_aux_size(num_bit_size, power, single_digit_bit_size);
    }

    proc pow_by_squaring(Slice* to, Slice* aux, CSlice num, umax power, Optims const& optims) -> Slice
    {
        assert(is_striped_number(num));
        //no two can alias
        assert(are_aliasing(*to, num) == false);
        assert(are_aliasing(*aux, num) == false);
        assert(are_aliasing(*aux, *to) == false);

        size_t bit_size = log2(num);
        const size_t required_size = required_pow_out_size(bit_size, power, DIGIT_BIT_SIZE);
        const size_t required_sigle_aux = required_pow_by_squaring_single_aux_swap_size(bit_size, power, DIGIT_BIT_SIZE);

        assert(to->size >= required_size);
        assert(aux->size >= required_pow_by_squaring_aux_size(bit_size, power, DIGIT_BIT_SIZE));

        if(power == 0 || (num.size == 1 && num[0] == 1))
        {
            Slice out = trim(*to, 1);
            out[0] = 1;
            return out;
        }

        if(num.size == 0)
            return trim(*to, 0);

        if(optims.mul_consts)
        {
            if(power == 1)
            {
                Slice trimmed_to = slice(*to, num.size);
                copy_slice(&trimmed_to, num);
                return trim(*to, num.size);
            }
            if(power == 2)
                return mul(to, aux, num, num, optims);
        }

        Slice output_aux1 = *to;
        Slice output_aux2 = trim(*aux, required_size);
        
        Slice after_out_aux = slice(*aux, required_size);

        Slice square_aux1 = slice_size(after_out_aux, required_sigle_aux*0, required_sigle_aux);
        Slice square_aux2 = slice_size(after_out_aux, required_sigle_aux*1, required_sigle_aux);

        Slice remianing_aux = slice(after_out_aux, 2*required_sigle_aux);

        //we switch storage on every asignment 
        // => we shoudl start with an appropiate one so that we finnish in the *to space
        //
        //lets say we do two assignments and start with aux1:
        // aux1 -> aux2 -> aux1    
        // => we should start with aux1 being *to 
        // => when we do even number of assignments we start with *to

        const size_t max_pos = find_last_set_bit(power);
        const size_t num_assignments = pop_count(power);
        if(num_assignments % 2 != 0)
            swap(&output_aux1, &output_aux2);

        CSlice curr_square = num;
        Slice curr_output = trim(output_aux1, 1);
        curr_output[0] = 1;

        for(size_t i = 0; i <= max_pos; i++)
        {
            
            Digit bit = (power >> i) & 1u;
            if(bit)
            {
                //curr_output *= curr_square;
                curr_output = mul(&output_aux2, &remianing_aux, curr_output, curr_square, optims);
                swap(&output_aux1, &output_aux2);
            }

            //if we are gonna stop next iteration dont waste time calculating the next mul
            if(i != max_pos)
            {
                //curr_square *= curr_square;
                curr_square = mul(&square_aux2, &remianing_aux, curr_square, curr_square, optims);
                swap(&square_aux1, &square_aux2);
            }
        }

        assert(curr_output.data == to->data);
        assert(is_striped_number(curr_output));
        return curr_output;
    }

    func required_pow_trivial_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size) -> size_t
    {
        if(power <= 1)
            return 0;

        //minus one since in the algorhitm we always multiply from output buffer to auxiliary
        // so that the final mutliply will be to the output buffer
        // => the maximum power that will be stored in the auxiliary is one less
        return required_pow_out_size(num_bit_size, power - 1, single_digit_bit_size);
    }
    
    func required_pow_aux_size(size_t num_bit_size, umax power, size_t single_digit_bit_size) -> size_t
    {
        return required_pow_trivial_aux_size(num_bit_size, power, single_digit_bit_size);
    }

    proc pow_trivial(Slice* to, Slice* aux, CSlice num, umax power, Optims const& optims) -> Slice
    {
        assert(is_striped_number(num));

        assert(are_aliasing(*to, num) == false);
        assert(are_aliasing(*aux, num) == false);
        assert(are_aliasing(*aux, *to) == false);

        const size_t required_to_size = required_pow_out_size(log2(num), power, DIGIT_BIT_SIZE);
        const size_t required_aux_size = required_pow_trivial_aux_size(log2(num), power, DIGIT_BIT_SIZE);
        assert(to->size >= required_to_size);
        assert(aux->size >= required_aux_size);

        if(power == 0 || (num.size == 1 && num[0] == 1))
        {
            Slice out = trim(*to, 1);
            out[0] = 1;
            return out;
        }

        if(num.size == 0)
            return trim(*to, 0);

        Slice output_aux1 = trim(*to, required_to_size);
        Slice output_aux2 = trim(*aux, required_aux_size);

        Slice remianing_aux = slice(*aux, required_aux_size);

        // we shoudl start with an appropiate one so that we finnish in the *to buffer
        const size_t num_assignments = power - 1;
        if(num_assignments % 2 != 0)
            swap(&output_aux1, &output_aux2);

        Slice curr = trim(output_aux1, num.size);
        copy_slice(&curr, num);

        for(size_t i = 0; i < num_assignments; i++)
        {
            Slice next = mul(&output_aux2, &remianing_aux, curr, num, optims);
            swap(&output_aux1, &output_aux2);
            curr = next;
        }

        assert(curr.data == to->data && "we should finish in the output buffer");
        assert(is_striped_number(curr));
        return curr;
    }

    proc pow(Slice* to, Slice* aux, CSlice num, umax power, Optims const& optims) -> Slice
    {
        size_t bitsize = log2(num);
        size_t at_least_to = required_pow_out_size(bitsize, power, DIGIT_BIT_SIZE);
        size_t at_least_aux = required_pow_trivial_aux_size(bitsize, power, DIGIT_BIT_SIZE);
        assert(to->size >= at_least_to);
        assert(aux->size >= at_least_aux);

        if(power < optims.pow_trivial_below_power)
            return pow_trivial(to, aux, num, power, optims);

        const size_t required_aux_size = required_pow_by_squaring_aux_size(bitsize, power, DIGIT_BIT_SIZE);

        if(aux->size >= required_aux_size)
            return pow_by_squaring(to, aux, num, power, optims);

        return pow_trivial(to, aux, num, power, optims);
    }

    func root_estimate_bit_size(size_t num_bit_size, size_t root) -> size_t {
        //const u64 upper_initial_estimate = cast(u64) 1 << ((log2x + n) / n);
        size_t estimate_bit_size = (num_bit_size + root) / root;
        return estimate_bit_size;
    }

    func required_root_out_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size) -> size_t {
        if(root == 0)
            return 1;

        if(root == 1)
            return div_round_up(num_bit_size, single_digit_bit_size);

        //the +1 is for the addition of the quotient we perform in place
        return div_round_up(root_estimate_bit_size(num_bit_size, root), single_digit_bit_size) + 1;
    }

    func optimal_root_aux_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size) -> size_t {

        size_t other_estimate_size = required_root_out_size(num_bit_size, root, single_digit_bit_size);
        size_t pow_aux_size = required_pow_by_squaring_aux_size(num_bit_size, root, single_digit_bit_size);
        size_t pow_to_size = required_pow_by_squaring_aux_size(num_bit_size, root, single_digit_bit_size);

        return other_estimate_size + pow_aux_size + pow_to_size;
    }

    func required_root_aux_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size) -> size_t {

        size_t other_estimate_size = required_root_out_size(num_bit_size, root, single_digit_bit_size);
        size_t pow_aux_size = required_pow_aux_size(num_bit_size, root, single_digit_bit_size);
        size_t pow_to_size = required_pow_aux_size(num_bit_size, root, single_digit_bit_size);

        return other_estimate_size + pow_aux_size + pow_to_size;
    }

    proc root(Slice* to, Slice* aux, CSlice num, size_t root, Optims const& optims) -> Slice
    {
        size_t num_bit_size = log2(num);   
        size_t required_to_size = required_root_out_size(num.size, root, DIGIT_BIT_SIZE);
        size_t required_aux_size = required_root_aux_size(num_bit_size, root, DIGIT_BIT_SIZE);

        assert(to->size >= required_to_size);
        assert(aux->size >= required_aux_size);

        assert(high_bits(cast(Digit) root) == 0 && "only small roots allowed - would require generalization of this algorhitm");
        if(num.size == 0)
        {
            to->data[0] = 1;
            return trim(*to, 1);
        }

        if(root <= 0)
        {
            Slice ret = trim(*to, num.size);
            copy_slice(&ret, num);
            return ret;
        }

        if(num.size == 1)
        {
            to->data[0] = cast(Digit) single_root_shifting(num[0], cast(Digit) root);
            return trim(*to, 1);
        }
        size_t estimate_bit_size = root_estimate_bit_size(num_bit_size, root);

        Slice out1 = *to;
        Slice out2 = trim(*aux, required_to_size);
        Slice rest_aux = slice(*aux, required_to_size);

        Slice initial_estimate = out2;

        null_slice(&initial_estimate);
        set_nth_bit(&initial_estimate, estimate_bit_size, 1);
        initial_estimate = striped_trailing_zeros(initial_estimate);
        //assert(is_striped_number(curr_estimate));

        Slice curr_estimate = initial_estimate;
        Slice prev_estimate = {};


        while(true)
        {
            prev_estimate = curr_estimate;
            assert(is_striped_number(prev_estimate));

            //new_upper_r = (root-1) * prev_estimate + num / single_pow_by_squaring(prev_estimate, root-1);

            
            //the sizes we will need to store the results:
            size_t powed_size = required_pow_out_size(log2(prev_estimate), root - 1, DIGIT_BIT_SIZE);
            size_t muled_size = required_mul_out_size(prev_estimate.size, 1);

            //assert(powed_size <= num.size && "is it?");
            
            //we will perform the operations in order of most memory consumption.
            // that way they can use as much still-not-occupied storage to speed up their execution
            Slice pow_to = trim(rest_aux, powed_size);
            Slice after_pow_aux = slice(rest_aux, powed_size);
            Slice powed = pow(&pow_to, &after_pow_aux, prev_estimate, root - 1, optims);

            size_t div_size = required_div_quotient_size(num.size, powed.size) + 1;
            size_t rem_size = required_div_remainder_size(num.size, powed.size);

            Slice div_to = trim(after_pow_aux, div_size);
            Slice after_div_aux = slice(after_pow_aux, div_size);
            Slice rem_to = trim(after_div_aux, rem_size);


            Div_Result res = div(&div_to, &rem_to, num, powed, optims);
            Slice dived = res.quotient;
            assert(is_striped_number(dived));

            size_t add_size = required_add_out_size(muled_size, dived.size);

            Slice new_estimate_to = out1;
            assert(new_estimate_to.data != prev_estimate.data && "storages must differ");
            assert(new_estimate_to.size >= add_size || true &&  
                "should fit into the prepared slot - maybe will fail because add_size is theoretical reuquired size"
                "but new_estimate_to.size is the actual real size which might - and probably be - smaller than the theoretical size"
                "we might have to lie about size to the algorhitm here!");

            //we perform the multiply add fused
            Batch_Overflow fused_result;
            
            if(dived.size != 0)
                fused_result = batch_fused_mul_add_overflow(&new_estimate_to, dived, prev_estimate, cast(Digit) root - 1, optims);
            else
                fused_result = batch_mul_overflow(&new_estimate_to, prev_estimate, cast(Digit) root - 1, optims);
                
            assert(fused_result.overflow == 0 && "should not overflow");
            assert(is_striped_number(fused_result.slice));
            swap(&out1, &out2);


            //r = new_upper_r / n;
            Slice new_upper_estimate = fused_result.slice;
            Slice new_estimate_pre = batch_div_overflow(&new_upper_estimate, new_upper_estimate, cast(Digit) root, optims).slice;
            Slice new_estimate = striped_trailing_zeros(new_estimate_pre);
            curr_estimate = new_estimate;

            //if(r >= prev_r)
                //break; //and the result is in prev_estimate (estimates shoudl be equal but just to be safe)
            int compared = compare(new_estimate, prev_estimate);
            if(compared >= 0)
                break;
        }

        if(prev_estimate.data == to->data)
            return prev_estimate;

        //if the current output is in the wrong storage copy over
        Slice dest = trim(*to, prev_estimate.size);
        copy_slice(&dest, prev_estimate);
        return dest;
    }

    struct Power_And_Powed
    {
        umax power = 0;
        umax powed = 0;

        bool constexpr operator ==(Power_And_Powed const&) const noexcept = default;
    };

    func highest_power_of_in_bits(size_t bit_size, size_t base) -> Power_And_Powed
    {
        assert(bit_size <= 64 && "no higher sizes supported");

        if(bit_size <= 1) 
            return {0, 1};

        switch(base)
        {
        case 10: 
            if(bit_size < 4)    return {0, 1};
            if(bit_size == 4)   return {1, 10};
            if(bit_size == 8)   return {2, 100};
            if(bit_size == 16)  return {4, 10'000};
            if(bit_size == 32)  return {9, 1'000'000'000};
            if(bit_size == 64)  return {19, 10'000'000'000'000'000'000};

            break;

        case 2:
            //if(bit_size == 64)
                //bit_size = 63;

            return {(bit_size - 1), 1ull << (bit_size - 1)};

        case 16:
            return {(bit_size - 1)/4, 1ull << (((bit_size - 1)/4) * 4)};
        }

        umax max = 0;
        if(bit_size < UMAX_BIT_SIZE)
            max = (1ull << bit_size) - 1;
        else
            max = cast(umax) -1;

        umax power = single_log_bsearch(max, base);
        umax powed = single_pow_by_squaring(base, power);

        return {power, powed};
    }

    func required_to_base_out_size(size_t num_size, size_t to_base, size_t single_digit_bit_size) -> size_t
    {
        Power_And_Powed highest = highest_power_of_in_bits(single_digit_bit_size / 2, to_base); 
        
        //into a single digit fits `power` powers of to_base
        // => we need to `power + 1` to_base to represent any number of the digit
        // => into the whole number fit `num_size * power` powers
        // => we need `num_size * (power + 1)` digits
        //@TODO: why do we need the +highest.power?
        return cast(size_t) (num_size * 2 * (highest.power + 1) + highest.power);
    }

    auto required_from_base_out_size(size_t rep_size, size_t from_base, size_t single_digit_bit_size) -> size_t
    {
        Power_And_Powed highest = highest_power_of_in_bits(single_digit_bit_size / 2, from_base); 

        //[][][][][]
        //----|----| 
        //  => 3 rep digits per num digit
        //  => we just need num digits / rep digits per digit rounded up
        return div_round_up(rep_size, cast(size_t) highest.power);
    }

    struct Optional_Slice
    {
        bool ok;
        Slice slice;
    };

    struct Pop_Digit
    {
        Slice buffer;
        CSlice number; //can also point into a different buffer from buffer
        size_t stashed_base = 0;
        size_t stashed_value = 0;
        size_t stashed_power = 0;
        Power_And_Powed highest;

        enum Error
        {
            OK,
            BASE_TOO_LARGE,
            BASE_TOO_SMALL,
            SMALL_BUFFER,
            NON_STRIPPED_NUM
        } error = OK;
        bool ok = true;
    };

    func pop_digit_init(Slice* buffer, Slice num, Digit base) -> Pop_Digit
    {
        Pop_Digit popper = Pop_Digit{*buffer, num, cast(size_t) base};

        popper.highest = highest_power_of_in_bits(DIGIT_HALF_BIT_SIZE, cast(size_t) base);
        if(high_bits(base) != 0)
        {
            popper.error = Pop_Digit::BASE_TOO_LARGE;
            popper.ok = false;
        }

        if(base < 2)
        {
            popper.error = Pop_Digit::BASE_TOO_SMALL;
            popper.ok = false;
        }

        if(buffer->size < num.size)
        {
            popper.error = Pop_Digit::SMALL_BUFFER;
            popper.ok = false;
        }

        if(is_striped_number(num) == false)
        {
            popper.error = Pop_Digit::NON_STRIPPED_NUM;
            popper.ok = false;
        }

        return popper;
    }

    proc pop_digit(Pop_Digit* state, Optims const& optims) -> Digit
    {
        if(state->stashed_power == 0)
        {
            if(state->number.size == 0)
            {
                state->ok = false; 
                return -1;
            }

            Batch_Overflow res = batch_div_overflow(&state->buffer, state->number, cast(Digit) state->highest.powed, optims);

            state->stashed_power = state->highest.power;
            state->stashed_value = res.overflow;
            state->number = striped_trailing_zeros(res.slice);
        }

        bool is_last_block = state->number.size == 0;
        if(is_last_block && state->stashed_value == 0)
        {
            state->ok = false; 
            return -1;
        }

        Digit digit = cast(Digit) state->stashed_value % state->stashed_base;
        state->stashed_value /= state->stashed_base;
        state->stashed_power -= 1;

        return digit;
    }

    struct Push_Digit
    {
        Slice buffer;
        size_t num_size = 0;
        umax stashed_base = 0;
        umax stashed_digit = 0;
        umax stashed_power = 0;
        Power_And_Powed highest;
    };

    func push_digit_init(Slice* buffer, size_t num_size = 0) -> Push_Digit
    {
        Push_Digit pusher = Push_Digit{*buffer, num_size};
        return pusher;
    }

    func push_digit(Push_Digit* state, Digit digit, Digit base, Optims const& optims, bool flush = false) -> bool
    {
        if(base < 2 || digit > base || state->num_size + 1 > state->buffer.size)
            return false;
        
        if(base != state->stashed_base)
        {
            if(state->stashed_power != 0)
                if(push_digit(state, 0, cast(Digit) state->stashed_base, optims, true))
                    return false;

            state->highest = highest_power_of_in_bits(DIGIT_BIT_SIZE, base);
            state->stashed_base = cast(size_t) base;
        }

        if(state->num_size + 1 > state->buffer.size)
            assert(false && "shouldnt be possible");

        //flush only attempts to flush stashed_digit into the num buffer
        if(flush == false)
        {
            state->stashed_digit *= base;
            state->stashed_digit += cast(umax) digit;
            state->stashed_power += 1;

            if(state->stashed_power < state->highest.power)
                return state;
        }

        Digit powed = cast(Digit) state->highest.powed;
        if(flush)
            powed = cast(Digit) single_pow_trivial(base, state->stashed_power);

        Slice curr_num = trim(state->buffer, state->num_size);
        size_t overflown_times = 0;

        Batch_Overflow mul_res = batch_mul_overflow(&curr_num, curr_num, powed, optims);
        if(mul_res.overflow != 0)
        {
            overflown_times += 1;
            state->num_size += 1;
            curr_num = trim(state->buffer, state->num_size);
            curr_num[state->num_size - 1] = mul_res.overflow;
        }

        Batch_Overflow add_res = batch_add_overflow_short(&curr_num, curr_num, cast(Digit) state->stashed_digit, Location::IN_PLACE);
        if(add_res.overflow != 0)
        {
            overflown_times += 1;
            state->num_size += 1;
            curr_num = trim(state->buffer, state->num_size);
            curr_num[state->num_size - 1] = add_res.overflow;
        }

        assert(overflown_times < 2 && "shouldt be possible to overflow 2 times");
        state->stashed_power = 0;
        state->stashed_digit = 0;
        return true;
    }

    proc push_digit_result(Push_Digit* state, Optims const& optims) -> Optional_Slice
    {
        if(state->stashed_power != 0)
        {
            if(push_digit(state, 0, cast(Digit) state->stashed_base, optims, true) == false)
                return Optional_Slice{false};
        }

        return Optional_Slice{true, trim(state->buffer, state->num_size)};
    }
    
    char to_char(Digit value)
    {
        constexpr char val_to_char_table[36] = {
            '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
            'U', 'V', 'W', 'X', 'Y', 'Z'
        };

        if(value > 36)
            return -1;
        else
            return val_to_char_table[value];
    }

    Digit to_digit(char value)
    {
        switch(value)
        {
            default: return -1;
            case '0': return 0;
            case '1': return 1;
            case '2': return 2;
            case '3': return 3;
            case '4': return 4;
            case '5': return 5;
            case '6': return 6;
            case '7': return 7;
            case '8': return 8;
            case '9': return 9;
            case 'A': return 10;
            case 'a': return 10;
            case 'B': return 11;
            case 'b': return 11;
            case 'C': return 12;
            case 'c': return 12;
            case 'D': return 13;
            case 'd': return 13;
            case 'E': return 14;
            case 'e': return 14;
            case 'F': return 15;
            case 'f': return 15;
            case 'G': return 16;
            case 'g': return 16;
            case 'H': return 17;
            case 'h': return 17;
            case 'I': return 18;
            case 'i': return 18;
            case 'J': return 19;
            case 'j': return 19;
            case 'K': return 20;
            case 'k': return 20;
            case 'L': return 21;
            case 'l': return 21;
            case 'M': return 22;
            case 'm': return 22;
            case 'N': return 23;
            case 'n': return 23;
            case 'O': return 24;
            case 'o': return 24;
            case 'P': return 25;
            case 'p': return 25;
            case 'Q': return 26;
            case 'q': return 26;
            case 'R': return 27;
            case 'r': return 27;
            case 'S': return 28;
            case 's': return 28;
            case 'T': return 29;
            case 't': return 29;
            case 'U': return 30;
            case 'u': return 30;
            case 'V': return 31;
            case 'v': return 32;
            case 'W': return 32;
            case 'w': return 32;
            case 'X': return 33;
            case 'x': return 33;
            case 'Y': return 34;
            case 'y': return 34;
            case 'Z': return 35;
            case 'z': return 35;
        }
    }

    proc from_chars(Slice* num, helpers::Slice<const char> rep, Digit base, Optims const& optims) -> Optional_Slice
    {
        assert(required_from_base_out_size(rep.size, base, DIGIT_BIT_SIZE) <= num->size 
            && "there must be enough size in num to represent the number even in the worst case scenario");

        Push_Digit state = push_digit_init(num);
        for(size_t i = 0; i < rep.size; i++)
        {
            Digit digit = to_digit(rep[i]);
            if(digit == -1)
                return {false};
                
            if(push_digit(&state, digit, base, optims) == false)
                return {false};
        }

        return push_digit_result(&state, optims);
    }

    //very fast to chars dumping
    proc to_chars(helpers::Slice<char>* rep, Slice* temp, CSlice num, Digit base, Optims const& optims) -> helpers::Slice<char>
    {
        assert(is_striped_number(num));
        assert(temp->size >= num.size && "temp must be big enough");
        assert(base >= 2 && "base must be bigger than two");
        assert(high_bits(base) == 0 && "base must be low bits only of Num type (so that short div algorhimt can be used)");

        size_t at_least = required_to_base_out_size(num.size, base, DIGIT_BIT_SIZE);
        assert(at_least <= rep->size && "there must be enough size rep represent the number even in the worst case scenario");

        bool is_first = true;
        size_t to_size = 0;

        if(num.size == 0)
            return trim(*rep, 0);

        Power_And_Powed highest = highest_power_of_in_bits(DIGIT_HALF_BIT_SIZE, base);
        assert(highest.powed > 1);

        CSlice val_current = num;
        bool iterate = true;
        while(iterate)
        {
            Batch_Overflow res = batch_div_overflow(temp, val_current, cast(Digit) highest.powed, optims);

            Digit highest_power_rem = res.overflow;
            val_current = striped_trailing_zeros(res.slice);

            bool is_last_block = val_current.size == 0;
            for(size_t j = 0; j < highest.power; j++)
            {
                if(is_last_block && highest_power_rem == 0)
                {
                    iterate = false;
                    break;
                }

                Digit rem = highest_power_rem % base;
                char digit = to_char(rem);
                (*rep)[to_size] = digit;
                to_size++;

                highest_power_rem /= base;
            }


            is_first = false;
        }

        helpers::Slice<char> converted = trim(*rep, to_size);

        size_t half_size = to_size / 2;
        for(size_t i = 0; i < half_size; i++)
            swap(&converted[i], &converted[to_size - i - 1]);

        return converted;
    }

    func single_factorial_calculate(size_t of_value) -> size_t
    {
        size_t calculated = 1;
        for(size_t factor = 2; factor <= of_value; factor++)
            calculated *= factor;

        return calculated;
    }

    static constexpr size_t MAX_NATIVE_FACTORIAL = 20;

    func single_factorial(size_t of_value) -> size_t
    {
        switch(of_value) {
            case 0: return single_factorial_calculate(0);
            case 1: return single_factorial_calculate(1);
            case 2: return single_factorial_calculate(2);
            case 3: return single_factorial_calculate(3);
            case 4: return single_factorial_calculate(4);
            case 5: return single_factorial_calculate(5);
            case 6: return single_factorial_calculate(6);
            case 7: return single_factorial_calculate(7);
            case 8: return single_factorial_calculate(8);
            case 9: return single_factorial_calculate(9);
            case 10: return single_factorial_calculate(10);
            case 11: return single_factorial_calculate(11);
            case 12: return single_factorial_calculate(12);
            case 13: return single_factorial_calculate(13);
            case 14: return single_factorial_calculate(14);
            case 15: return single_factorial_calculate(15);
            case 16: return single_factorial_calculate(16);
            case 17: return single_factorial_calculate(17);
            case 18: return single_factorial_calculate(18);
            case 19: return single_factorial_calculate(19);
            case 20: return single_factorial_calculate(20);
            default: return -1;
        }
    }

    func required_factorial_size(size_t of_value) -> size_t
    {
        if(of_value <= MAX_NATIVE_FACTORIAL)
            return 1;
        else
            return of_value - MAX_NATIVE_FACTORIAL;
    }

    proc factorial(Slice* out, size_t of_value, Optims const& optims) -> Slice
    {
        assert(out->size >= required_factorial_size(of_value));

        bool did_finish = true;
        if(of_value < MAX_NATIVE_FACTORIAL)
        {
            size_t fact = single_factorial(of_value);
            return from_number(out, fact);
        }

        const size_t highest_possible_native = single_factorial(MAX_NATIVE_FACTORIAL);
        Slice current_value = from_number(out, highest_possible_native);
        for(size_t factor = MAX_NATIVE_FACTORIAL + 1; factor <= of_value; factor++)
            current_value = mul(out, current_value, cast(Digit) factor, optims);

        return current_value;
    }

    #undef proc
    #undef func
    #undef cast
}

#undef tiny
#endif // TINY_BIG_NUM_ONCE