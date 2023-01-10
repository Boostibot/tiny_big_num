#include "tiny_big_num.h"

#define cast(...) (__VA_ARGS__)
//#define CHAR_BIT 8

//#ifdef INCLUDED_NEW_TYPE
//#ifndef TINY_NUM_LIB_NAME
//#ifdef INCLUDED_TINY_NUM_TYPE_DEF
//#define TINY_NUM_LIB_NAME tiny_num
//#else
//#error "name must be provided!"
//#endif
//#endif 
//
//#define CONCAT_(a, b) a ## b
//#define CONCAT(a, b) CONCAT_(a, b)
//
////name macro prepends the lib name to every exported symbol
//#define N(a) CONCAT(TINY_NUM_LIB_NAME, a)
//#endif

Optims make_def_optims()
{
    Optims optims;
    optims.mul_shift = DO_OPTIM_MUL_SHIFT;
    optims.mul_half_bits = DO_OPTIM_MUL_HALF_BITS;
    optims.div_shift = DO_OPTIM_DIV_SHIFT;

    optims.max_recursion_depth = MAX_RECURSION_DEPTH;
    optims.mul_quadratic_both_below_size = MUL_QUADRATIC_BOTH_BELOW_SIZE;
    optims.mul_quadratic_single_below_size = MUL_QUADRATIC_SINGLE_BELOW_SIZE;
    optims.pow_trivial_below_power = TRIVIAL_POW_BELOW_POWER;

    return optims;
}

enum Iter_Direction 
{
    FORWARD,
    BACKWARD
};

typedef struct Single_Overflow
{
    Digit value;
    Digit overflow;
} Single_Overflow;

typedef struct Batch_Overflow
{
    Slice slice;
    Digit overflow;
} Batch_Overflow;

typedef struct Power_And_Powed
{
    umax power;
    umax powed;
} Power_And_Powed;

Result ok_result(Slice output)
{
    Result res = {OK, output};
    return res;
}

Result error_result(State state)
{
    Result res = {state};
    return res;
}

Div_Result ok_div_result(Slice quot, Slice rem)
{
    Div_Result res = {OK, quot, rem};
    return res;
}

Div_Result error_div_result(State state)
{
    Div_Result res = {state};
    return res;
}

Char_Result error_char_result(State state)
{
    Char_Result res = {state};
    return res;
}

inline size_t max(size_t a, size_t b)
{
    return a > b ? a : b;
}

inline size_t min(size_t a, size_t b)
{
    return a < b ? a : b;
}

inline size_t div_round_up(size_t value, size_t to_multiple_of)
{
    return (value + to_multiple_of - 1) / to_multiple_of;
}

size_t pop_count_u32(uint32_t val)
{
    uint32_t i = cast(uint32_t) val;
    i = i - ((i >> 1) & 0x55555555);                    // add pairs of bits
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);     // quads
    i = (i + (i >> 4)) & 0x0F0F0F0F;                    // groups of 8
    return cast(size_t) (i * 0x01010101) >> 24;         // horizontal sum of bytes
}

size_t find_last_set_bit(umax val) 
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

size_t log2(umax val)
{
    return find_last_set_bit(val);
}

size_t find_first_set_bit(umax val)
{
    if(val == 0)
        return 0;

    //10110000
    //&   1111=> +4 and shift
    //00001011
    //&     11=> +0
    //00001011
    //&      1=> +0
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

size_t pop_count(umax val)
{
    uint32_t low = cast(uint32_t) (val & 0xFFFF'FFFF);
    uint32_t high = cast(uint32_t) (val >> 32);

    return pop_count_u32(low) + pop_count_u32(high);
}

umax single_pow_trivial(umax x, umax n)
{
    umax res = 1;
    for(umax i = 0; i < n; i++)
        res *= x;

    return res;
}

umax single_pow_by_squaring(umax x, umax n)
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

umax single_root_shifting(umax x, umax n)
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

umax single_root_newton(umax of_value, umax power)
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

umax single_log_heyley(umax of_value, umax base)
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

umax single_log_bsearch(umax of_value, umax base)
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

CSlice to_cslice(Slice items)
{
    CSlice out = {items.data, items.size};
    return out;
}

//asserts the index is in range and returns it
// this compiles to noop on release mode and is used because
// we have to hide more statemnets into expressions acessed through macros
// see at( , ) macro
inline size_t pass_through_bounds_check(size_t index, size_t bounds)
{
    assert(index < bounds && "index out of bounds!");
    return index;
}

#define at(items, i) \
    (items).data[pass_through_bounds_check(i, (items).size)]

Slice slice_range(Slice items, size_t from, size_t to)
{
    assert(to >= from && "must be a valid range");
    assert(from <= items.size && to <= items.size && "sliced portion must be entirely within the given range");
    Slice out = {items.data + from, to - from};
    return out;
}

Slice slice_size(Slice items, size_t from, size_t count) {
    return slice_range(items, from, from + count);
}

Slice slice(Slice items, size_t from) {
    return slice_range(items, from, items.size);
}

Slice trim(Slice items, size_t to_size) {
    return slice_range(items, 0, to_size);
}

CSlice cslice_range(CSlice items, size_t from, size_t to)
{
    assert(to >= from && "must be a valid range");
    assert(from <= items.size && to <= items.size && "sliced portion must be entirely within the given range");
    CSlice out = {items.data + from, to - from};
    return out;
}

CSlice cslice_size(CSlice items, size_t from, size_t count) {
    return cslice_range(items, from, from + count);
}

CSlice cslice(CSlice items, size_t from) {
    return cslice_range(items, from, items.size);
}

CSlice ctrim(CSlice items, size_t to_size) {
    return cslice_range(items, 0, to_size);
}

Digit last(CSlice slice)
{
    assert(slice.size > 0 && "must not be empty!");
    return slice.data[slice.size - 1];
}

#define swap(a, b, T)   \
    {                   \
        T temp = *a;    \
        *a = *b;        \
        *b = temp;      \
    }                   \

void swap_slice(Slice* a, Slice* b)      { swap(a, b, Slice); }
void swap_cslice(CSlice* a, CSlice* b)   { swap(a, b, CSlice); }
void swap_size(size_t* a, size_t* b)     { swap(a, b, size_t); }
void swap_digit(Digit* a, Digit* b)      { swap(a, b, Digit); }

#undef swap

bool are_aliasing(CSlice left, CSlice right)
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

bool are_one_way_aliasing(CSlice before, CSlice after)
{ 
    return (before.data + before.size > after.data) && (after.data > before.data);
}

//into cpp
void copy_slice(Slice* to, CSlice from)
{
    assert(to->size == from.size && "sizes must match");
    memmove(to->data, from.data, from.size * sizeof(Digit));
}

void null_slice(Slice* to)
{
    memset(to->data, 0, to->size * sizeof(Digit));
}

static const size_t CHAR_BIT = 8; //@TEMP

static const size_t DIGIT_BIT_SIZE = sizeof(Digit) * CHAR_BIT;
static const size_t UMAX_BIT_SIZE = sizeof(umax) * CHAR_BIT;
static const size_t DIGIT_HALF_BIT_SIZE = DIGIT_BIT_SIZE / 2;
static const size_t MAX_NATIVE_FACTORIAL = 20; //@TODO

Digit set_bit(Digit field, size_t bit_pos, Digit val) {
    assert(val == 0 || val == 1);
    assert(bit_pos < DIGIT_BIT_SIZE);

    const Digit bit = val << bit_pos;
    const Digit mask = cast(Digit) 1 << bit_pos;
    const Digit res = (field & ~mask) | bit;
    return res;
};

Digit get_bit(Digit field, size_t bit_pos) 
{
    return (field >> bit_pos) & 1u;
};

Digit set_nth_bit(Slice* num, size_t i, Digit val)
{
    const size_t digit_i = i / DIGIT_BIT_SIZE;
    const size_t bit_i = i % DIGIT_BIT_SIZE;

    assert(digit_i < num->size);

    Digit* digit = &at(*num, digit_i);
    *digit = set_bit(*digit, bit_i, val);
};

Digit get_nth_bit(CSlice num, size_t i) {
    const size_t digit_i = i / DIGIT_BIT_SIZE;
    const size_t bit_i = i % DIGIT_BIT_SIZE;

    return get_bit(at(num, digit_i), bit_i);
};

Digit high_mask_i(size_t index) {
    assert(index < DIGIT_BIT_SIZE);
    const Digit full_mask = -1;
    return cast(Digit) (full_mask << index);
}

Digit low_mask_i(size_t index) {
    return cast(Digit) ~high_mask_i(index);
}

Digit high_bits_i(Digit value, size_t index) {
    assert(index < DIGIT_BIT_SIZE);
    return cast(Digit) (value >> index);
}

Digit low_bits_i(Digit value, size_t index) {
    return cast(Digit) (value & low_mask_i(index));
};  

Digit combine_bits_i(Digit low, Digit high, size_t index) {
    assert(index < DIGIT_BIT_SIZE);
    return cast(Digit) (low_bits_i(low, index) | (high << index));
};

Digit dirty_combine_bits_i(Digit low, Digit high, size_t index) {
    assert(index < DIGIT_BIT_SIZE);
    assert(high_bits_i(low, index) == 0 && "low must not have high bits use combine_bits instead");
    return cast(Digit) (low | (high << index));
};


Digit high_mask()                               { return high_mask_i(DIGIT_HALF_BIT_SIZE); }
Digit low_mask()                                { return low_mask_i(DIGIT_HALF_BIT_SIZE); }
Digit high_bits(Digit value)                    { return high_bits_i(value, DIGIT_HALF_BIT_SIZE); }
Digit low_bits(Digit value)                     { return low_bits_i(value, DIGIT_HALF_BIT_SIZE); }
Digit combine_bits(Digit low, Digit high)       { return combine_bits_i(low, high, DIGIT_HALF_BIT_SIZE); }
Digit dirty_combine_bits(Digit low, Digit high) { return dirty_combine_bits_i(low, high, DIGIT_HALF_BIT_SIZE); }

size_t find_last_set_digit(CSlice num)
{
    size_t i = num.size;
    for(; i-- > 0;)
        if(at(num, i) != 0)
            break;

    return i;
}

size_t find_first_set_digit(CSlice num)
{
    size_t i = 0;
    for(; i < num.size; i++)
        if(at(num, i) != 0)
            break;

    return i;
}

Slice striped_trailing_zeros(Slice num)
{
    return trim(num, find_last_set_digit(to_cslice(num)) + 1);
}


CSlice cstriped_trailing_zeros(CSlice num)
{
    return ctrim(num, find_last_set_digit(num) + 1);
}

bool is_striped_number(CSlice num)
{
    if(num.size == 0)
        return true;

    return num.data[num.size - 1] != 0;
}

umax to_number(CSlice bignum) 
{
    size_t digit_i = find_last_set_digit(bignum);
    CSlice stripped = ctrim(bignum, digit_i);

    const size_t total_size = stripped.size * sizeof(Digit);
    assert(total_size <= sizeof(umax));

    umax out = 0;
    for(size_t i = 0; i < stripped.size; i++)
    {
        size_t shift_by = i * DIGIT_BIT_SIZE;
        assert(shift_by < UMAX_BIT_SIZE);
        Digit curr = at(stripped, i);
        out |= cast(umax) curr << shift_by;
    }

    return out;
}

Slice from_number(Slice* bignum, umax from)
{
    constexpr size_t count = sizeof(umax) / sizeof(Digit);
    assert(bignum->size >= count);

    size_t i = 0;
    for(; i < count; i++)
    {
        umax curr_digit = from >> (i * DIGIT_BIT_SIZE);
        if(curr_digit == 0)
            break;

        at(*bignum, i) = cast(Digit) curr_digit;
    }

    return trim(*bignum, i);
}



Digit single_add_no_overflow(Digit a1, Digit a2, Digit a3)
{
    #ifdef NDEBUG
    return a1 + a2 + a3;
    #else
    const Digit arr[] = {a1, a2, a3};
    Digit sum = a1;

    for(size_t i = 1; i < 3; i++)
    {
        Digit new_sum = sum + arr[i];
        assert(new_sum >= sum && "must not overflow");
        sum = new_sum;
    }

    return sum;
    #endif // NDEBUG
}

Single_Overflow single_add_overflow_any(Digit carry, Digit a1, Digit a2, Digit a3)
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

Single_Overflow single_sub_overflow_any(Digit carry, Digit left, Digit a2, Digit a3)
{
    const Digit res_low = high_mask() + low_bits(left) - (low_bits(a2) + low_bits(a3)) - carry;
    const Digit middle_carry = low_mask() - high_bits(res_low);

    const Digit res_high = high_mask() + high_bits(left) - (high_bits(a2) + high_bits(a3)) - middle_carry;
    const Digit new_carry = low_mask() - high_bits(res_high);

    return Single_Overflow{combine_bits(res_low, res_high), new_carry};
}

Single_Overflow single_add_overflow(Digit left, Digit right, Digit carry)
{
    return single_add_overflow_any(carry, left, right, 0);
}

Single_Overflow single_sub_overflow(Digit left, Digit right, Digit carry)
{
    return single_sub_overflow_any(carry, left, right, 0);
}

Single_Overflow single_complement_overflow(Digit val, Digit carry)
{
    const Digit val_inv = ~val;

    const Digit complement_low = low_bits(val_inv) + carry;
    const Digit middle_carry = high_bits(complement_low);

    const Digit complement_high = high_bits(val_inv) + middle_carry;
    const Digit new_carry = high_bits(complement_high);

    return Single_Overflow{combine_bits(complement_low, complement_high), new_carry};
}

Single_Overflow single_shift_up_overflow(Digit low_item, size_t by_bits, Digit high_item)
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

    const Digit low = high_bits_i(low_item, remaining_bits);
    const Digit high = low_bits_i(high_item, remaining_bits);

    const Digit composed = low | high << by_bits;

    const Digit out = high_bits_i(high_item, remaining_bits);
    const Digit shifted_out = out;
    return Single_Overflow{composed, shifted_out};
}

Single_Overflow single_shift_down_overflow(Digit low_item, size_t by_bits, Digit high_item)
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


    const Digit low = high_bits_i(low_item, by_bits);
    const Digit high = low_bits_i(high_item, by_bits);

    const Digit composed = low | high << remaining_bits;

    const Digit out = low_bits_i(low_item, by_bits);
    const Digit shifted_out = out << remaining_bits;
    return Single_Overflow{composed, shifted_out};
}

enum Mul_Overflow_Optims
{
    NONE,
    HIGH_BITS_ONLY,
    LOW_BITS_ONLY,
};

Single_Overflow single_mul_overflow(Digit left, Digit right, Digit last_value, Mul_Overflow_Optims mul_optims)
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

Single_Overflow single_div_overflow(Digit left, Digit right, Digit carry_in)
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

bool single_is_power_of_two(Digit num)
{
    return num != 0 && (num & (num - 1)) == 0;
}



Batch_Overflow batch_add_or_sub_overflow_short(Slice* to, CSlice left, Digit carry, bool is_addition, size_t from)
{
    assert(to->size >= left.size);
    assert(are_one_way_aliasing(left, to_cslice(*to)) == false);

    const Slice trimmed_to = trim(*to, left.size);
    size_t j = from;
    for (; carry != 0 && j < left.size != 0; j++)
    {
        const Digit digit = at(left, j);
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

        at(trimmed_to, j) = patch_res;
    }

    if(to->data == left.data)
    {
        Slice remainign_to = slice(trimmed_to, j);
        copy_slice(&remainign_to, cslice(left, j));
    }

    return Batch_Overflow{trimmed_to, carry};
}

Batch_Overflow batch_add_overflow_short(Slice* to, CSlice left, Digit right, size_t from)
{
    return batch_add_or_sub_overflow_short(to, left, right, true, from);
}

Batch_Overflow batch_sub_overflow_short(Slice* to, CSlice left, Digit right, size_t from)
{
    return batch_add_or_sub_overflow_short(to, left, right, false, from);
}

Batch_Overflow batch_add_overflow_long(Slice* to, CSlice left, CSlice right, Digit carry_in)
{
    assert(to->size >= right.size);
    assert(to->size >= left.size);
    assert(are_one_way_aliasing(left, to_cslice(*to)) == false);
    assert(high_bits(carry_in) == 0);

    if(left.size < right.size)
        swap_cslice(&left, &right);

    Digit carry = carry_in;
    for (size_t i = 0; i < right.size; i++)
    {
        Single_Overflow res = single_add_overflow(at(left, i), at(right, i), carry);
        at(*to, i) = res.value;
        carry = res.overflow;
    }

    return batch_add_overflow_short(to, left, carry, right.size);
}

Batch_Overflow batch_sub_overflow_long(Slice* to, CSlice left, CSlice right, Digit carry_in)
{
    assert(to->size >= right.size);
    assert(to->size >= left.size);
    assert(are_one_way_aliasing(left, to_cslice(*to)) == false);
    assert(carry_in == 0 || carry_in == 1);

    Digit carry = carry_in;
    const size_t min_size = min(left.size, right.size);
    const size_t max_size = max(left.size, right.size);
    for (size_t i = 0; i < min_size; i++)
    {
        Single_Overflow oveflow = single_sub_overflow(at(left, i), at(right, i), carry);
        at(*to, i) = oveflow.value;
        carry = oveflow.overflow;
    }

    for (size_t i = left.size; i < right.size; i++)
    {
        Single_Overflow oveflow = single_sub_overflow(cast(Digit) 0, at(right, i), carry);
        at(*to, i) = oveflow.value;
        carry = oveflow.overflow;
    }

    if(right.size < left.size)
        return batch_sub_overflow_short(to, left, carry, right.size);

    return Batch_Overflow{trim(*to, max_size), carry};
}

Batch_Overflow batch_complement_overflow(Slice* to, CSlice left, Digit carry_in)
{
    assert(to->size >= left.size);
    assert(is_striped_number(left));
    assert(are_one_way_aliasing(left, to_cslice(*to)) == false);
    assert(carry_in == 1 || carry_in == 0);

    Digit carry = carry_in;
    for(size_t i = 0; i < left.size; i++)
    {
        Single_Overflow res = single_complement_overflow(at(left, i), carry);
        at(*to, i) = res.value;
        carry = res.overflow;
    }

    return Batch_Overflow{trim(*to, left.size), carry};
}

//We allow shifting while iterating in both directions:
// 
// FOR SHIFT_UP: the FORWARD direction is useful in general case
//  (for example when shifting array in place)
// and the BACKWARD direction can be used while shifting by more then one element
//  and writes to already passed memory are neccessary
// 
// FOR SHIFT_DOWN: its the polar opposite
Batch_Overflow batch_shift_up_overflow(Slice* out, CSlice in, size_t by_bits, 
    Iter_Direction direction, Digit carry_in)
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
        assert(are_one_way_aliasing(in, to_cslice(*out)) == false);
        Digit prev = carry_in;
        for (size_t i = 0; i < in.size; i++)
        {
            Digit curr = at(in, i);
            Single_Overflow res = single_shift_up_overflow(prev, by_bits, curr);
            at(*out, i) = res.value;
            prev = curr;
        }

        Single_Overflow carry_out = single_shift_up_overflow(prev, by_bits, 0);
        return Batch_Overflow{out_trimmed, carry_out.value};
    }
    else
    {
        assert(are_one_way_aliasing(to_cslice(*out), in) == false);
        const Digit shifted_out = last(in);
        Digit prev = shifted_out;
        for (size_t i = in.size; i-- > 1; )
        {
            Digit curr = at(in, i - 1);
            at(*out, i) = single_shift_up_overflow(curr, by_bits, prev).value;
            prev = curr;
        }

        {
            Single_Overflow res = single_shift_up_overflow(carry_in, by_bits, prev);
            at(*out, 0) = res.value;
        }


        Single_Overflow carry_out = single_shift_up_overflow(shifted_out, by_bits, 0);
        return Batch_Overflow{out_trimmed, carry_out.value};
    }
}

Batch_Overflow batch_shift_down_overflow(Slice* out, CSlice in, size_t by_bits, 
    Iter_Direction direction, Digit carry_in)
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
        assert(are_one_way_aliasing(in, to_cslice(*out)) == false);

        const Digit shifted_out = at(in, 0);
        Digit prev = shifted_out;
        for (size_t i = 0; i < in.size - 1; i++)
        {
            Digit curr = at(in, i + 1);
            Single_Overflow res = single_shift_down_overflow(prev, by_bits, curr);
            at(*out, i) = res.value;
            prev = curr;
        }

        {
            Single_Overflow res = single_shift_down_overflow(prev, by_bits, carry_in);
            at(*out, in.size - 1) = res.value;
        }

        Single_Overflow carry_out = single_shift_down_overflow(0, by_bits, shifted_out);
        return Batch_Overflow{out_trimmed, carry_out.value};
    }
    else
    {
        assert(are_one_way_aliasing(to_cslice(*out), in) == false);

        Digit prev = carry_in;
        for (size_t i = in.size; i-- > 0;)
        {
            Digit curr = at(in, i);
            at(*out, i) = single_shift_down_overflow(curr, by_bits, prev).value;
            prev = curr;
        }

        Single_Overflow carry_out = single_shift_down_overflow(0, by_bits, prev);
        return Batch_Overflow{out_trimmed, carry_out.value};
    }
}

Batch_Overflow batch_mul_overflow(Slice* to, CSlice left, Digit right, const Optims* optims, Digit carry_in)
{
    assert(to->size >= left.size);
    assert(are_one_way_aliasing(left, to_cslice(*to)) == false);

    Digit carry = carry_in;
    Slice trimmed_to = trim(*to, left.size);

    if(right == 0) 
        return Batch_Overflow{trim(*to, 0), 0};

    if(right == 1)
    {
        if(trimmed_to.data != left.data)
            copy_slice(&trimmed_to, left);
        return Batch_Overflow{trimmed_to, 0};
    }

    if(optims->mul_shift)
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

    #define do_carry_loop(operation)                    \
        {                                               \
            for (size_t i = 0; i < left.size; i++)      \
            {                                           \
                const Single_Overflow res = operation;  \
                at(*to, i) = res.value;                 \
                carry = res.overflow;                   \
            }                                           \
        }                                               \

    if(optims->mul_half_bits && high_bits(right) == 0)
        do_carry_loop((single_mul_overflow(at(left, i), right, carry, Mul_Overflow_Optims::LOW_BITS_ONLY)))
    else if(optims->mul_half_bits && low_bits(right) == 0)
        do_carry_loop((single_mul_overflow(at(left, i), right, carry, Mul_Overflow_Optims::HIGH_BITS_ONLY)))
    else
        do_carry_loop((single_mul_overflow(at(left, i), right, carry, Mul_Overflow_Optims::NONE)))

        #undef do_carry_loop
        return Batch_Overflow{trimmed_to, carry};
}

Batch_Overflow batch_div_overflow(Slice* to, CSlice left, Digit right, const Optims* optims, Digit carry_in)
{
    assert(to->size >= left.size);
    assert(high_bits(right) == 0 && "only works for divisors under half bit size");
    assert(right != 0 && "cannot divide by zero");
    assert(are_one_way_aliasing(to_cslice(*to), left) == false);

    Slice trimmed_to = trim(*to, left.size);

    if(right == 1)
    {
        if(trimmed_to.data != left.data)
            copy_slice(&trimmed_to, left);
        return Batch_Overflow{trimmed_to, 0}; 
    }

    if(optims->div_shift)
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
        Single_Overflow res = single_div_overflow(at(left, i), right, carry);
        at(*to, i) = res.value;
        carry = res.overflow;
    }

    return Batch_Overflow{trimmed_to, carry};
}

Digit batch_rem_overflow(CSlice left, Digit right, const Optims* optims, Digit carry_in)
{
    assert(high_bits(right) == 0 && "only works for divisors under half bit size");
    assert(right != 0 && "cannot divide by zero");

    if(single_is_power_of_two(right))
    {
        size_t shift_bits = find_last_set_bit(right);
        if(left.size == 0)
            return carry_in;

        return low_bits_i(at(left, 0), shift_bits);
    }

    Digit carry = carry_in;
    for (size_t i = left.size; i-- > 0; )
    {
        Single_Overflow res = single_div_overflow(at(left, i), right, carry);
        carry = res.overflow;
    }

    return carry;
}


//return -1 if left < right
//        0 if left == right
//        1 if left > right   
int compare(CSlice left, CSlice right)
{
    assert(is_striped_number(left));
    assert(is_striped_number(right));
    if(left.size < right.size)
        return -1;

    if(left.size > right.size)
        return 1;

    for(size_t i = left.size; i-- > 0;)
    {
        Digit l = at(left, i);
        Digit r = at(right, i);
        if(l < r)
            return -1;
        if(l > r)
            return 1;
    }

    return 0;
}

bool is_equal(CSlice left, CSlice right)
{
    if(left.size != right.size)
        return false;

    for(size_t i = 0; i < left.size; i++)
        if(at(left, i) != at(right, i))
            return false;

    return true;
}

size_t required_add_out_size(size_t left_size, size_t right_size)
{
    return max(left_size, right_size) + 1;
}

Result add_(Slice* to, CSlice left, CSlice right)
{
    assert(is_striped_number(left));
    assert(is_striped_number(right));
    if(to->size < required_add_out_size(left.size, right.size))
        return error_result(OUT_SIZE_TOO_SMALL);

    Batch_Overflow res = batch_add_overflow_long(to, left, right);
    if(res.overflow == 0)
    {
        assert(is_striped_number(to_cslice(res.slice)));
        return ok_result(res.slice);
    }

    size_t slice_size = res.slice.size;
    Slice trimmed_to = trim(*to, slice_size + 1);
    at(trimmed_to, slice_size) = res.overflow;

    assert(is_striped_number(to_cslice(trimmed_to)));

    return ok_result(trimmed_to);
}

size_t required_sub_out_size(size_t left_size, size_t right_size)
{
    return left_size;
}

Result sub_(Slice* to, CSlice left, CSlice right)
{
    assert(is_striped_number(left));
    assert(is_striped_number(right));
    if(to->size < required_sub_out_size(left.size, right.size))
        return error_result(OUT_SIZE_TOO_SMALL);

    Batch_Overflow res = batch_sub_overflow_long(to, left, right);
    
    if(res.overflow != 0)
    {
        Result result = {OK_EXEPTIONAL, res.slice};
        return result;
    }

    Slice output = striped_trailing_zeros(res.slice);
    return ok_result(output);
}

size_t required_mul_out_size(size_t left_size, size_t right_size)
{
    return left_size + right_size;
}

Result mul_short(Slice* to, CSlice left, Digit right, const Optims* optims)
{
    assert(is_striped_number(left));
    if(to->size < required_mul_out_size(left.size, 1))
        return error_result(OUT_SIZE_TOO_SMALL);

    Batch_Overflow res = batch_mul_overflow(to, left, right, optims, 0);
    if(res.overflow == 0)
    {
        assert(is_striped_number(to_cslice(res.slice)));
        return ok_result(res.slice);
    }

    size_t slice_size = res.slice.size;
    Slice trimmed_to = trim(*to, slice_size + 1);
    at(trimmed_to, slice_size) = res.overflow;

    assert(is_striped_number(to_cslice(trimmed_to)));
    return ok_result(trimmed_to);
}

size_t required_div_quotient_size(size_t num_size, size_t den_size)
{
    if(num_size < den_size)
        return 0;

    return num_size - den_size + 1;
}

size_t required_div_remainder_size(size_t num_size, size_t den_size)
{
    return min(den_size + 1, num_size);;
}

Div_Short_Result div_short(Slice* to, CSlice left, Digit right, const Optims* optims)
{
    assert(is_striped_number(left));
    if(right == 0)
    {
        Div_Short_Result result = {UNDEFINED_VALUE};
        return result;
    }

    if(to->size < required_mul_out_size(left.size, 1))
    {
        Div_Short_Result result = {OUT_SIZE_TOO_SMALL};
        return result;
    }

    Batch_Overflow batch = batch_div_overflow(to, left, right, optims, 0);
    Slice quotient = striped_trailing_zeros(batch.slice);
    Div_Short_Result result = {OK, quotient, batch.overflow};
    return result;
}

Div_Short_Result div_short_in_place(Slice* left, Digit right, const Optims* optims)
{
    return div_short(left, to_cslice(*left), right, optims);
}

Div_Result div_bit_by_bit_(Slice* quotient, Slice* remainder, CSlice num, CSlice den, const Optims* optims, bool DO_QUOTIENT)
{
    assert(is_striped_number(num));
    assert(is_striped_number(den));
    assert(are_aliasing(to_cslice(*quotient), num) == false);
    assert(are_aliasing(to_cslice(*remainder), num) == false);
    assert(are_aliasing(to_cslice(*quotient), den) == false);
    assert(are_aliasing(to_cslice(*remainder), den) == false);

    const size_t required_quot_size = DO_QUOTIENT ? required_div_quotient_size(num.size, den.size) : 0;
    const size_t required_rem_size = required_div_remainder_size(num.size, den.size);

    if(DO_QUOTIENT)
        assert(are_aliasing(to_cslice(*quotient), to_cslice(*remainder)) == false);

    assert(den.size > 0 && "cannot divide by zero!");

    if(num.size == 0)
    {
        Slice stripped_quotient = trim(*quotient, 0);
        Slice stripped_remainder = trim(*remainder, 0);
        return ok_div_result(stripped_quotient, stripped_remainder);
    }

    if(num.size < den.size)
    {
        assert(remainder->size >= num.size);
        Slice stripped_remainder = trim(*remainder, num.size);
        Slice stripped_quotient = trim(*quotient, 0);

        copy_slice(&stripped_remainder, num);
        return ok_div_result(stripped_quotient, stripped_remainder);
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
            at(trimmed_remainder, 0) = remainder_val;
            remainder_size = 1;
        }

        Slice stripped_remainder = trim(trimmed_remainder, remainder_size);
        return ok_div_result(stripped_quotient, stripped_remainder);
    }

    Slice curr_remainder = trim(trimmed_remainder, 0);
    for(size_t i = DIGIT_BIT_SIZE * num.size; i-- > 0; ) 
    {
        size_t last_size = curr_remainder.size;
        Batch_Overflow shift_res = batch_shift_up_overflow(&curr_remainder, to_cslice(curr_remainder), 1);
        const Digit num_ith_bit = get_nth_bit(num, i);

        if(shift_res.overflow != 0 || (last_size == 0 && num_ith_bit == 1))
        {
            curr_remainder = trim(trimmed_remainder, last_size + 1);
            at(curr_remainder, last_size) = shift_res.overflow;
        }

        set_nth_bit(&trimmed_remainder, 0, num_ith_bit);
        if(compare(to_cslice(curr_remainder), den) >= 0)
        {
            Batch_Overflow sub_res = batch_sub_overflow_long(&curr_remainder, to_cslice(curr_remainder), den);
            curr_remainder = striped_trailing_zeros(curr_remainder);
            assert(sub_res.overflow == 0 && "should not overflow");

            if(DO_QUOTIENT)
                set_nth_bit(&trimmed_quotient, i, 1);
        }
    }

    assert(is_striped_number(to_cslice(curr_remainder)));
    Slice stripped_quotient = DO_QUOTIENT 
        ? striped_trailing_zeros(trimmed_quotient)
        : trim(trimmed_quotient, 0);

    return ok_div_result(stripped_quotient, curr_remainder);
}


Div_Result div_bit_by_bit(Slice* quotient, Slice* remainder, CSlice num, CSlice den, const Optims* optims)
{
    return div_bit_by_bit_(quotient, remainder, num, den, optims, true);
}

Div_Result div(Slice* quotient, Slice* remainder, CSlice num, CSlice den, const Optims* optims)
{
    return div_bit_by_bit(quotient, remainder, num, den, optims);
}

Result rem(Slice* remainder, CSlice num, CSlice den, const Optims* optims)
{
    Slice quotient = {nullptr, 0};
    Div_Result div_res = div_bit_by_bit_(&quotient, remainder, num, den, optims, false);
    Result result = {div_res.state, div_res.remainder};

    assert(div_res.quotient.size == 0 && "in case of only remainder the quotient size should be 0");
    assert(div_res.quotient.data == nullptr && "and data nullptr as we set it");
    return result;
}

Batch_Overflow batch_fused_mul_add_overflow(Slice* to, CSlice added, CSlice multiplied, Digit coeficient, 
    const Optims* optims, Digit add_carry, Digit mul_carry)
{   
    assert(added.size >= multiplied.size);
    assert(to->size >= added.size);
    assert(are_one_way_aliasing(added, to_cslice(*to)) == false);
    assert(are_one_way_aliasing(multiplied, to_cslice(*to)) == false);

    Slice trimmed_to = trim(*to, added.size);
    if(coeficient == 0) 
    {
        if (to->data == added.data)
            return Batch_Overflow{trim(*to, 0), 0};
        else
        {
            copy_slice(&trimmed_to, added);
            return Batch_Overflow{trimmed_to, 0};
        }
    }

    if(coeficient == 1)
        return batch_add_overflow_long(to, added, multiplied);

    Slice trimmed_to = trim(*to, added.size);
    if(optims->mul_shift 
        && single_is_power_of_two(coeficient))
    {
        size_t shift_bits = find_last_set_bit(coeficient);
        const Slice trimmed_to = trim(*to, added.size);
        Digit prev_muled = mul_carry;
        for (size_t j = 0; j < multiplied.size; j++)
        {
            const Digit curr_muled = at(multiplied, j);
            const Digit curr_added = at(added, j);

            Single_Overflow mul_res = single_shift_up_overflow(prev_muled, shift_bits, curr_muled);
            Single_Overflow add_res = single_add_overflow(curr_added, mul_res.value, add_carry);

            at(trimmed_to, j) = add_res.value;

            prev_muled = curr_muled;
            add_carry = add_res.overflow;
        }

        mul_carry = single_shift_up_overflow(prev_muled, shift_bits, 0).value;
    }
    else
    {
        for (size_t j = 0; j < multiplied.size; j++)
        {
            const Digit curr_muled = at(multiplied, j);
            const Digit curr_added = at(added, j);

            Single_Overflow mul_res = single_mul_overflow(curr_muled, coeficient, mul_carry);
            Single_Overflow add_res = single_add_overflow(curr_added, mul_res.value, add_carry);

            at(trimmed_to, j) = add_res.value;

            mul_carry = mul_res.overflow;
            add_carry = add_res.overflow;
        }
    }

    const Digit combined_carry = single_add_no_overflow(mul_carry, add_carry, 0);
    return batch_add_overflow_short(&trimmed_to, added, combined_carry, multiplied.size);
}

size_t required_mul_quadratic_out_size(size_t left_size, size_t right_size)
{
    return required_mul_out_size(left_size, right_size);
}

size_t required_mul_quadratic_aux_size(size_t left_size, size_t right_size)
{
    return required_mul_out_size(left_size, right_size);
}

Result mul_quadratic(Slice* to, Slice* temp, CSlice left, CSlice right, const Optims* optims)
{
    assert(are_aliasing(to_cslice(*to), left) == false);
    assert(are_aliasing(to_cslice(*to), right) == false);
    assert(are_aliasing(to_cslice(*temp), left) == false);
    assert(are_aliasing(to_cslice(*temp), right) == false);
    assert(are_aliasing(to_cslice(*to), to_cslice(*temp)) == false);

    if(left.size < right.size)
        swap_cslice(&left, &right);

    const size_t max_result_size = required_mul_out_size(left.size, right.size);
    const size_t max_temp_size = required_mul_out_size(left.size, right.size);

    if(to->size < max_result_size)
        return error_result(OUT_SIZE_TOO_SMALL);
    if(temp->size < max_temp_size)
        return error_result(AUX_SIZE_TOO_SMALL);

    Slice trimmed_to = trim(*to, max_result_size);
    Slice trimmed_temp = trim(*temp, max_temp_size);
    null_slice(&trimmed_to);

    for(size_t i = 0; i < right.size; i++)
    {
        Batch_Overflow mul_res = batch_mul_overflow(&trimmed_temp, left, at(right, i), optims);

        //batch_mul_overflow is limited to left.size elems => in case of overflow
        // even if the overflow would fit into temp it is not added => we add it back
        size_t mul_size = mul_res.slice.size;
        Slice mul_slice = trim(trimmed_temp, mul_size + 1); 
        at(mul_slice, mul_size) = mul_res.overflow;

        Slice shifted_to = slice(trimmed_to, i);
        assert(shifted_to.size >= mul_slice.size);

        Batch_Overflow add_res = batch_add_overflow_long(&shifted_to, to_cslice(shifted_to), to_cslice(mul_slice));
        assert(add_res.overflow == 0 && "no final carry should be left");
    }

    Slice output = striped_trailing_zeros(trimmed_to);
    return ok_result(output);
}

size_t required_mul_quadratic_fused_out_size(size_t left_size, size_t right_size)
{
    return required_mul_out_size(left_size, right_size);
}

Result mul_quadratic_fused(Slice* to, CSlice left, CSlice right, const Optims* optims)
{
    const size_t max_result_size = required_mul_out_size(left.size, right.size);
    assert(are_aliasing(to_cslice(*to), left) == false);
    assert(are_aliasing(to_cslice(*to), right) == false);
    if(to->size < max_result_size)
        return error_result(OUT_SIZE_TOO_SMALL);

    if(left.size < right.size)
        swap_cslice(&left, &right);


    Slice trimmed_to = trim(*to, max_result_size);

    //we null only the place for the first multiplication / addition
    // the rest is nulled implicitly through the overflow handling code
    // this saves us one additinal iteration through the data
    Slice first_iter = trim(*to, left.size);
    null_slice(&first_iter);

    for(size_t i = 0; i < right.size; i++)
    {
        Slice shifted_to = slice_size(trimmed_to, i, left.size);

        Batch_Overflow mul_add_res = batch_fused_mul_add_overflow(&shifted_to, to_cslice(shifted_to), left, at(right, i), optims, 0, 0);
        at(trimmed_to, i + left.size) = mul_add_res.overflow;
    }

    Slice output = striped_trailing_zeros(trimmed_to);
    return ok_result(output);
}

size_t required_mul_karatsuba_out_size(size_t left_size, size_t right_size)
{
    return required_mul_out_size(left_size, right_size);
}

size_t required_mul_karatsuba_aux_size(size_t left_size, size_t right_size, size_t recursion_depth)
{
    if(left_size < right_size)
        swap_size(&left_size, &right_size);

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

Result mul_(Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims, size_t depth);

Result mul_karatsuba_(Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims, size_t depth, bool is_run_alone)
{
    assert(is_striped_number(left));
    assert(is_striped_number(right));

    assert(are_aliasing(to_cslice(*to), left) == false);
    assert(are_aliasing(to_cslice(*aux), left) == false);
    assert(are_aliasing(to_cslice(*to), right) == false);
    assert(are_aliasing(to_cslice(*aux), right) == false);
    assert(are_aliasing(to_cslice(*aux), to_cslice(*to)) == false);

    CSlice x = left;
    CSlice y = right;

    //Perfomrs the karatsuba multiplicaton step
    //x*y = (x1*b + x2)(y1*b + y2) = (x1*y1)b^2 + (x1*y2 + x2*y1)b + (x2*y2)
    //    = z1*b^2 + z2*b + z3

    // z1 = x1*y1
    // z3 = x2*y2
    // z2 = x1*y2 + x2*y1 = (x1 + x2)(y1 + y2) - z1 - z3

    if(x.size < y.size)
        swap_cslice(&x, &y);

    if(is_run_alone)
    {
        if(x.size == 0 || y.size == 0)
            return ok_result(trim(*to, 0));

        if(y.size == 1)
            return mul_short(to, x, at(y, 0), optims);
    }

    //size_t split_digits = min(x.size / 2, y.size / 2);
    //we split according to the bigger one. This might seem odd but in the case where we will split too much
    // and y1 will be 0 we will have one less multiplication to worry about
    size_t base = max(x.size, y.size) / 2; 
    size_t capped_digits = min(y.size, base);
    assert(base != 0);

    CSlice x1 = cslice(x, base);
    CSlice x2 = cstriped_trailing_zeros(ctrim(x, base));

    CSlice y1 = cslice(y, capped_digits);
    CSlice y2 = cstriped_trailing_zeros(ctrim(y, capped_digits));


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

    assert(are_aliasing(to_cslice(z1_slot), to_cslice(z3_slot)) == false);

    //we multyply into z3 and null the area between it and begginign of z1 
    // (ie the place where z2 will get added to)

    Slice z3 = mul_(&z3_slot, &remaining_aux, x2, y2, optims, depth + 1).output;
    Slice z3_to_z1 = slice_range(trimmed_out, z3.size, z1_from_index);
    null_slice(&z3_to_z1);

    //we multiply into z1 and null the remaining size to the end of trimmed_out
    // this way trimmed out should look like the following:
    // 000[ z1 ]000[ z3 ]
    Slice z1 = mul_(&z1_slot, &remaining_aux, x1, y1, optims, depth + 1).output;
    Slice z1_up = slice(trimmed_out, z1.size + z1_from_index);
    null_slice(&z1_up);

    size_t used_to = 0;
    Slice z2_slot = slice_size(*aux, used_to, required_z2_size);
    used_to += required_z2_size;

    Slice x_sum_slot = slice_size(*aux, used_to, required_add_x_sum_size);
    used_to += required_add_x_sum_size;

    Slice y_sum_slot = slice_size(*aux, used_to, required_add_y_sum_size);
    used_to += required_add_y_sum_size;

    Slice x_sum = add(&x_sum_slot, x1, x2).output;
    Slice y_sum = add(&y_sum_slot, y1, y2).output;

    Slice x_sum_y_sum = mul_(&z2_slot, &remaining_aux, to_cslice(x_sum), to_cslice(y_sum), optims, depth + 1).output;
    Slice x_sum_y_sum_m_z1 = sub_(&x_sum_y_sum, to_cslice(x_sum_y_sum), to_cslice(z1)).output;
    Slice z2 = sub_(&x_sum_y_sum_m_z1, to_cslice(x_sum_y_sum_m_z1), to_cslice(z3)).output;

    //instead of multiplying z2 by base we add it to the appropriate position
    Slice out_z2_up = slice(trimmed_out, base);

    Batch_Overflow add_res = batch_add_overflow_long(&out_z2_up, to_cslice(out_z2_up), to_cslice(z2));
    assert(add_res.overflow == 0 && "should not overflow");

    Slice output = striped_trailing_zeros(trimmed_out);
    return ok_result(output);
}

Result mul_karatsuba(Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims)
{
    return mul_karatsuba_(to, aux, left, right, optims, 0, true);
}

size_t optimal_mul_aux_size(size_t left_size, size_t right_size, size_t recursion_depth)
{
    return required_mul_karatsuba_aux_size(left_size, right_size, recursion_depth);
}

Result mul_(Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims, size_t depth)
{
    assert(is_striped_number(left));
    assert(is_striped_number(right));
    assert(are_aliasing(to_cslice(*to), left) == false);
    assert(are_aliasing(to_cslice(*to), right) == false);

    if(to->size < required_mul_out_size(left.size, right.size))
        return error_result(OUT_SIZE_TOO_SMALL);

    if(left.size < right.size)
        swap_cslice(&left, &right);

    if(right.size == 0)
        return ok_result(trim(*to, 0));

    if(right.size == 1)
        return mul_short(to, left, at(right, 0), optims);

    if(right.size >= optims->mul_quadratic_single_below_size 
        && left.size >= optims->mul_quadratic_both_below_size 
        && depth < optims->max_recursion_depth)
    {
        size_t required_aux = required_mul_karatsuba_aux_size(left.size, right.size, 1);
        if(aux->size >= required_aux)
        {
            Result result = mul_karatsuba_(to, aux, left, right, optims, depth, false);
            assert(is_striped_number(to_cslice(result.output)));
            return result;
        }
    }

    Result result = mul_quadratic_fused(to, left, right, optims);
    assert(is_striped_number(to_cslice(result.output)));
    return result;
}

Result mul(Slice* to, Slice* aux, CSlice left, CSlice right, const Optims* optims)
{
    return mul_(to, aux, left, right, optims, 0);
}

size_t bit_size(CSlice left)
{
    assert(is_striped_number(left));
    if(left.size == 0)
        return 0;

    size_t last_log = find_last_set_bit(last(left));
    return last_log + (left.size - 1) * DIGIT_BIT_SIZE;
}

size_t log2(CSlice left)
{
    return bit_size(left);
}

size_t required_pow_out_size(size_t num_bit_size, umax power, size_t digit_bit_size)
{
    if(power == 0 || num_bit_size == 0)
        return 1;

    size_t bit_size = (num_bit_size + 1) * power;
    size_t item_size = div_round_up(bit_size, digit_bit_size) + 1;

    return item_size;
}

size_t required_pow_by_squaring_single_aux_swap_size(size_t num_bit_size, umax power, size_t digit_bit_size)
{
    if(power == 0)
        return 0;

    size_t iters = find_last_set_bit(power);
    size_t bits = (num_bit_size + 1);
    //we square the number on every iteration => the number of bits doubles => 2^iters == 1 << iters
    size_t bit_size = bits << iters; // == bits * 2^iters == bits * 2^[log2(power)] ~~ bits*power
    size_t item_size = div_round_up(bit_size, digit_bit_size) + 1;
    return item_size;
}

size_t required_pow_by_squaring_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size)
{
    if(power == 0)
        return 0;

    size_t square = required_pow_by_squaring_single_aux_swap_size(num_bit_size, power, digit_bit_size);
    size_t out_swap = required_pow_out_size(num_bit_size, power, digit_bit_size);

    return square * 2 + out_swap;
}

size_t optimal_pow_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size) {
    return required_pow_by_squaring_aux_size(num_bit_size, power, digit_bit_size);
}

Result pow_by_squaring(Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims)
{
    assert(is_striped_number(num));
    //no two can alias
    assert(are_aliasing(to_cslice(*to), num) == false);
    assert(are_aliasing(to_cslice(*aux), num) == false);
    assert(are_aliasing(to_cslice(*aux), to_cslice(*to)) == false);

    size_t bit_size = log2(num);
    const size_t required_size = required_pow_out_size(bit_size, power, DIGIT_BIT_SIZE);
    const size_t required_sigle_aux = required_pow_by_squaring_single_aux_swap_size(bit_size, power, DIGIT_BIT_SIZE);

    if(to->size < required_size)
        return error_result(OUT_SIZE_TOO_SMALL);

    if(aux->size < required_pow_by_squaring_aux_size(bit_size, power, DIGIT_BIT_SIZE))
        return error_result(AUX_SIZE_TOO_SMALL);

    if(power == 0 || (num.size == 1 && at(num, 0) == 1))
    {
        Slice out = trim(*to, 1);
        at(out, 0) = 1;
        return ok_result(out);
    }

    if(num.size == 0)
    {
        Slice out = trim(*to, 0);
        return ok_result(out);
    }

    //@REMOVE
    #if 0
    if(optims->mul_consts)
    {
        if(power == 1)
        {
            Slice trimmed_to = slice(*to, num.size);
            copy_slice(&trimmed_to, num);
            return ok_result(trimmed_to);
        }
        if(power == 2)
            return mul(to, aux, num, num, optims);
    }
    #endif

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
    // aux1
    // => we should start with aux1 being *to 
    // => when we do even number of assignments we start with *to

    const size_t max_pos = find_last_set_bit(power);
    const size_t num_assignments = pop_count(power);
    if(num_assignments % 2 != 0)
        swap_slice(&output_aux1, &output_aux2);

    CSlice curr_square = num;
    Slice curr_output = trim(output_aux1, 1);
    at(curr_output, 0) = 1;

    for(size_t i = 0; i <= max_pos; i++)
    {

        Digit bit = (power >> i) & 1u;
        if(bit)
        {
            //curr_output *= curr_square;
            Result mul_res = mul(&output_aux2, &remianing_aux, to_cslice(curr_output), curr_square, optims);
            curr_output = mul_res.output;
            assert(mul_res.state == OK);
            swap_slice(&output_aux1, &output_aux2);
        }

        //if we are gonna stop next iteration dont waste time calculating the next mul
        if(i != max_pos)
        {
            //curr_square *= curr_square;
            Result mul_res = mul(&square_aux2, &remianing_aux, curr_square, curr_square, optims);
            curr_square = to_cslice(mul_res.output);
            assert(mul_res.state == OK);
            swap_slice(&square_aux1, &square_aux2);
        }
    }

    assert(curr_output.data == to->data);
    assert(is_striped_number(to_cslice(curr_output)));
    return ok_result(curr_output);
}

size_t required_pow_trivial_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size)
{
    if(power <= 1)
        return 0;

    //minus one since in the algorhitm we always multiply from output buffer to auxiliary
    // so that the final mutliply will be to the output buffer
    // => the maximum power that will be stored in the auxiliary is one less
    return required_pow_out_size(num_bit_size, power - 1, digit_bit_size);
}

size_t required_pow_aux_size(size_t num_bit_size, umax power, size_t digit_bit_size)
{
    return required_pow_trivial_aux_size(num_bit_size, power, digit_bit_size);
}

Result pow_trivial(Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims)
{
    assert(is_striped_number(num));

    assert(are_aliasing(to_cslice(*to), num) == false);
    assert(are_aliasing(to_cslice(*aux), num) == false);
    assert(are_aliasing(to_cslice(*aux), to_cslice(*to)) == false);

    const size_t required_to_size = required_pow_out_size(log2(num), power, DIGIT_BIT_SIZE);
    const size_t required_aux_size = required_pow_trivial_aux_size(log2(num), power, DIGIT_BIT_SIZE);
    if(to->size >= required_to_size)
        return error_result(OUT_SIZE_TOO_SMALL);
    if(aux->size >= required_aux_size)
        return error_result(AUX_SIZE_TOO_SMALL);

    if(power == 0 || (num.size == 1 && at(num, 0) == 1))
    {
        Slice out = trim(*to, 1);
        at(out, 0) = 1;
        return ok_result(out);
    }

    if(num.size == 0)
    {
        Slice out = trim(*to, 0);
        return ok_result(out);
    }

    Slice output_aux1 = trim(*to, required_to_size);
    Slice output_aux2 = trim(*aux, required_aux_size);

    Slice remianing_aux = slice(*aux, required_aux_size);

    // we shoudl start with an appropiate one so that we finnish in the *to buffer
    const size_t num_assignments = power - 1;
    if(num_assignments % 2 != 0)
        swap_slice(&output_aux1, &output_aux2);

    Slice curr = trim(output_aux1, num.size);
    copy_slice(&curr, num);

    for(size_t i = 0; i < num_assignments; i++)
    {
        Result mul_res = mul(&output_aux2, &remianing_aux, to_cslice(curr), num, optims);
        Slice next = mul_res.output;
        assert(mul_res.state == OK);
        swap_slice(&output_aux1, &output_aux2);
        curr = next;
    }

    assert(curr.data == to->data && "we should finish in the output buffer");
    assert(is_striped_number(to_cslice(curr)));
    return ok_result(curr);
}

Result pow(Slice* to, Slice* aux, CSlice num, umax power, const Optims* optims)
{
    size_t bitsize = log2(num);
    size_t at_least_to = required_pow_out_size(bitsize, power, DIGIT_BIT_SIZE);
    size_t at_least_aux = required_pow_trivial_aux_size(bitsize, power, DIGIT_BIT_SIZE);
    if(to->size < at_least_to) 
        return error_result(OUT_SIZE_TOO_SMALL);
    if(aux->size < at_least_aux)
        return error_result(OUT_SIZE_TOO_SMALL);

    if(power < optims->pow_trivial_below_power)
        return pow_trivial(to, aux, num, power, optims);

    const size_t required_aux_size = required_pow_by_squaring_aux_size(bitsize, power, DIGIT_BIT_SIZE);

    if(aux->size >= required_aux_size)
        return pow_by_squaring(to, aux, num, power, optims);

    return pow_trivial(to, aux, num, power, optims);
}

size_t root_estimate_bit_size(size_t num_bit_size, umax root){
    //const u64 upper_initial_estimate = cast(u64) 1 << ((log2x + n) / n);
    umax estimate_bit_size = (num_bit_size + root) / root;
    return cast(size_t) estimate_bit_size;
}

size_t required_root_out_size(size_t num_bit_size, umax root, size_t digit_bit_size){
    if(root == 0)
        return 1;

    if(root == 1)
        return div_round_up(num_bit_size, digit_bit_size);

    //the +1 is for the addition of the quotient we perform in place
    return div_round_up(root_estimate_bit_size(num_bit_size, root), digit_bit_size) + 1;
}

size_t optimal_root_aux_size(size_t num_bit_size, umax root, size_t digit_bit_size){

    size_t other_estimate_size = required_root_out_size(num_bit_size, root, digit_bit_size);
    size_t pow_aux_size = optimal_pow_aux_size(num_bit_size, root, digit_bit_size);
    size_t pow_to_size = required_pow_out_size(num_bit_size, root, digit_bit_size);

    return other_estimate_size + pow_aux_size + pow_to_size;
}

size_t required_root_aux_size(size_t num_bit_size, umax root, size_t digit_bit_size){

    size_t other_estimate_size = required_root_out_size(num_bit_size, root, digit_bit_size);
    size_t pow_aux_size = required_pow_aux_size(num_bit_size, root, digit_bit_size);
    size_t pow_to_size = required_pow_out_size(num_bit_size, root, digit_bit_size);

    return other_estimate_size + pow_aux_size + pow_to_size;
}

Result root(Slice* to, Slice* aux, CSlice num, umax root, const Optims* optims)
{
    size_t num_bit_size = log2(num);   
    size_t required_to_size = required_root_out_size(num.size, root, DIGIT_BIT_SIZE);
    size_t required_aux_size = required_root_aux_size(num_bit_size, root, DIGIT_BIT_SIZE);

    if(to->size < required_to_size)
        return error_result(OUT_SIZE_TOO_SMALL);

    if(aux->size < required_aux_size)
        return error_result(AUX_SIZE_TOO_SMALL);

    assert(high_bits(cast(Digit) root) == 0 && "only small roots allowed - would require generalization of this algorhitm");
    if(num.size == 0)
    {
        Slice out = trim(*to, 0);
        return ok_result(out);
    }

    if(root <= 0)
        return error_result(UNDEFINED_VALUE);

    if(num.size == 1)
    {
        Slice out = trim(*to, 1);
        at(out, 0) = cast(Digit) single_root_shifting(at(num, 0), cast(Digit) root);
        return ok_result(out);
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
        assert(is_striped_number(to_cslice(prev_estimate)));

        //new_upper_r = (root-1) * prev_estimate + num / single_pow_by_squaring(prev_estimate, root-1);


        //the sizes we will need to store the results:
        size_t powed_size = required_pow_out_size(log2(to_cslice(prev_estimate)), root - 1, DIGIT_BIT_SIZE);
        size_t muled_size = required_mul_out_size(prev_estimate.size, 1);

        //assert(powed_size <= num.size && "is it?");

        //we will perform the operations in order of most memory consumption.
        // that way they can use as much still-not-occupied storage to speed up their execution
        Slice pow_to = trim(rest_aux, powed_size);
        Slice after_pow_aux = slice(rest_aux, powed_size);

        Result pow_result = pow(&pow_to, &after_pow_aux, to_cslice(prev_estimate), root - 1, optims);
        assert(pow_result.state == OK && "shouldnt fail");

        Slice powed = pow_result.output;
        size_t div_size = required_div_quotient_size(num.size, powed.size) + 1;
        size_t rem_size = required_div_remainder_size(num.size, powed.size);

        Slice div_to = trim(after_pow_aux, div_size);
        Slice after_div_aux = slice(after_pow_aux, div_size);
        Slice rem_to = trim(after_div_aux, rem_size);

        Div_Result res = div(&div_to, &rem_to, num, to_cslice(powed), optims);
        assert(res.state == OK && "shouldnt fail");

        Slice dived = res.quotient;
        assert(is_striped_number(to_cslice(dived)));

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
            fused_result = batch_fused_mul_add_overflow(&new_estimate_to, to_cslice(dived), to_cslice(prev_estimate), cast(Digit) root - 1, optims, 0, 0);
        else
            fused_result = batch_mul_overflow(&new_estimate_to, to_cslice(prev_estimate), cast(Digit) root - 1, optims);

        assert(fused_result.overflow == 0 && "should not overflow");
        assert(is_striped_number(to_cslice(fused_result.slice)));
        swap_slice(&out1, &out2);


        //r = new_upper_r / n;
        Slice new_upper_estimate = fused_result.slice;
        Slice new_estimate_pre = batch_div_overflow(&new_upper_estimate, to_cslice(new_upper_estimate), cast(Digit) root, optims).slice;
        Slice new_estimate = striped_trailing_zeros(new_estimate_pre);
        curr_estimate = new_estimate;

        //if(r >= prev_r)
        //break; //and the result is in prev_estimate (estimates shoudl be equal but just to be safe)
        int compared = compare(to_cslice(new_estimate), to_cslice(prev_estimate));
        if(compared >= 0)
            break;
    }

    if(prev_estimate.data == to->data)
        return ok_result(prev_estimate);

    //if the current output is in the wrong storage copy over
    Slice dest = trim(*to, prev_estimate.size);
    copy_slice(&dest, to_cslice(prev_estimate));
    return ok_result(dest);
}


Power_And_Powed highest_power_of_in_bits(size_t bit_size, umax base)
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

size_t required_to_base_out_size(size_t num_size, umax to_base, size_t digit_bit_size)
{
    Power_And_Powed highest = highest_power_of_in_bits(digit_bit_size / 2, to_base); 

    //into a single digit fits `power` powers of to_base
    // => we need to `power + 1` to_base to represent any number of the digit
    // => into the whole number fit `num_size * power` powers
    // => we need `num_size * (power + 1)` digits
    //@TODO: why do we need the +highest.power?
    return cast(size_t) (num_size * 2 * (highest.power + 1) + highest.power);
}

size_t required_from_base_out_size(size_t rep_size, umax from_base, size_t digit_bit_size)
{
    Power_And_Powed highest = highest_power_of_in_bits(digit_bit_size / 2, from_base); 

    //[][][][][]
    //----|----| 
    //  => 3 rep digits per num digit
    //  => we just need num digits / rep digits per digit rounded up
    return div_round_up(rep_size, cast(size_t) highest.power);
}

typedef struct Pop_Digit
{
    size_t num_size;
    umax stashed_base;
    umax stashed_value;
    umax stashed_power;
    Power_And_Powed highest;
    State state;
} Pop_Digit;

State pop_digit_can_continue(const Pop_Digit* pop)
{
    return pop->state;
}

Pop_Digit pop_digit_init(Slice* out, CSlice num, Digit base)
{
    assert(is_striped_number(num));

    Pop_Digit pop = {num.size, cast(umax) base};

    pop.highest = highest_power_of_in_bits(DIGIT_HALF_BIT_SIZE, cast(size_t) base);
    pop.state = OK;

    if(high_bits(base) != 0)
        pop.state = ARGUMENTS_TOO_BIG;
    else if(base < 2)
        pop.state = ARGUMENTS_TOO_SMALL;
    else if(num.size == 0)
        pop.state = OK_EXEPTIONAL;
    else if(out->size < num.size)
        pop.state = OUT_SIZE_TOO_SMALL;
    else if(num.data != out->data)
    {
        Slice used_portion = trim(*out, num.size);
        copy_slice(&used_portion, num);
    }

    return pop;
}

Digit pop_digit(Pop_Digit* pop, Slice* out, const Optims* optims)
{
    if(pop->state != OK)
        return -1;

    if(pop->stashed_power == 0)
    {
        if(pop->num_size == 0)
        {
            pop->state = OUT_SIZE_TOO_SMALL; 
            return -1;
        }

        Slice used_portion = trim(*out, pop->num_size);
        Batch_Overflow res = batch_div_overflow(out, to_cslice(used_portion), cast(Digit) pop->highest.powed, optims);

        pop->stashed_power = pop->highest.power;
        pop->stashed_value = res.overflow;
        pop->num_size = find_first_set_digit(to_cslice(res.slice)) + 1;
    }

    bool is_last_block = pop->num_size == 0;
    if(is_last_block && pop->stashed_value == 0)
    {
        pop->state = OUT_SIZE_TOO_SMALL; 
        return -1;
    }

    Digit digit = cast(Digit) pop->stashed_value % pop->stashed_base;
    pop->stashed_value /= pop->stashed_base;
    pop->stashed_power -= 1;

    if(pop->num_size == 0)
        pop->state = OK_EXEPTIONAL; 

    return digit;
}

Result pop_digit_result(Pop_Digit* pop, Slice* out, const Optims* optims)
{
    Slice trimmed_out = trim(*out, pop->num_size);
    if(pop->stashed_power == 0)
        return ok_result(trimmed_out);

    Result muled = mul_short(out, to_cslice(trimmed_out), cast(Digit) pop->stashed_power, optims);
    if(muled.state != OK)
        return muled;
    
    Result added = add_short(out, to_cslice(trimmed_out), cast(Digit) pop->stashed_value);
    return added;
}


typedef struct Push_Digit
{
    size_t num_size;
    umax stashed_base;
    umax stashed_digit;
    umax stashed_power;
    Power_And_Powed highest;
    State state;
} Push_Digit;

State push_digit_can_continue(const Push_Digit* push)
{
    return push->state;
}

Push_Digit push_digit_init(Slice* out, CSlice num, Digit base)
{
    assert(is_striped_number(num));
    Push_Digit push = {num.size, cast(umax) base};

    push.highest = highest_power_of_in_bits(DIGIT_BIT_SIZE, base);
    push.state = OK;

    if(base < 2)
        push.state = ARGUMENTS_TOO_SMALL;
    else if(out->size < num.size)
        push.state = OUT_SIZE_TOO_SMALL;
    else if(num.data != out->data && num.size != 0)
    {
        Slice used_portion = trim(*out, num.size);
        copy_slice(&used_portion, num);
    }

    return push;
}

State push_digit_or_flush(Push_Digit* push, Slice* out, Digit digit, const Optims* optims, bool flush)
{
    if(digit > push->stashed_base)
    {
        push->state = ARGUMENTS_TOO_BIG;
        return ARGUMENTS_TOO_BIG;
    }

    if(push->num_size + 1 > out->size)
        return OUT_SIZE_TOO_SMALL;

    //flush only attempts to flush stashed_digit into the num buffer
    if(flush == false)
    {
        push->stashed_digit *= push->stashed_base;
        push->stashed_digit += cast(umax) digit;
        push->stashed_power += 1;

        if(push->stashed_power < push->highest.power)
            return OK;
    }

    Digit powed = cast(Digit) push->highest.powed;
    if(flush)
        powed = cast(Digit) single_pow_trivial(push->stashed_base, push->stashed_power);

    Slice curr_num = trim(*out, push->num_size);
    size_t overflown_times = 0;

    Result mul_res = mul_short(out, to_cslice(curr_num), powed, optims);
    Result add_res = add_short(out, to_cslice(curr_num), cast(Digit) push->stashed_digit);

    size_t new_size = add_res.output.size;
    assert(new_size - push->num_size < 2 && "shouldt be possible to overflow 2 times");

    push->stashed_power = 0;
    push->stashed_digit = 0;
    push->num_size = new_size;
    return OK;
}

State push_digit(Push_Digit* push, Slice* out, Digit digit, const Optims* optims)
{
    return push_digit_or_flush(push, out, digit, optims, false);
}

Result push_digit_result(Push_Digit* push, Slice* out, const Optims* optims)
{
    if(push->stashed_power != 0)
    {
        State flush_state = push_digit_or_flush(push, out, 0, optims, true);
        if(flush_state != OK)
            return error_result(flush_state);
    }

    Slice out = trim(*out, push->num_size);
    return ok_result(*out);
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

Result from_chars(Slice* out, Char_CSlice rep, Digit base, const Optims* optims)
{
    if(required_from_base_out_size(rep.size, base, DIGIT_BIT_SIZE) <= out->size)
        return error_result(OUT_SIZE_TOO_SMALL);

    CSlice empty = {0};
    Push_Digit push = push_digit_init(out, empty, base);
    for(size_t i = 0; i < rep.size; i++)
    {
        Digit digit = to_digit(at(rep, i));
        if(digit == -1)
            return error_result(UNDEFINED_VALUE);

        State state = push_digit(&push, out, digit, optims);
        if(state != OK)
            return error_result(state);
    }

    return push_digit_result(&push, out, optims);
}

Char_Result to_chars(Char_Slice* rep, Slice* aux, CSlice num, Digit base, const Optims* optims)
{
    assert(is_striped_number(num));
    size_t at_least = required_to_base_out_size(num.size, base, DIGIT_BIT_SIZE);

    if(at_least <= rep->size)
        error_char_result(OUT_SIZE_TOO_SMALL);

    if(aux->size >= num.size)
        error_char_result(AUX_SIZE_TOO_SMALL);

    if(base <= 2)
        error_char_result(ARGUMENTS_TOO_SMALL);
        
    if(high_bits(base) != 0)    
        error_char_result(ARGUMENTS_TOO_BIG);
        

    bool is_first = true;
    size_t to_size = 0;

    if(num.size == 0)
    {
        Char_Slice out = {rep->data, 0};
        Char_Result result = {OK, out};
        return result;
    }

    Power_And_Powed highest = highest_power_of_in_bits(DIGIT_HALF_BIT_SIZE, base);
    assert(highest.powed > 1);

    CSlice val_current = num;
    bool iterate = true;
    while(iterate)
    {
        Batch_Overflow res = batch_div_overflow(aux, val_current, cast(Digit) highest.powed, optims);

        Digit highest_power_rem = res.overflow;
        val_current = to_cslice(striped_trailing_zeros(res.slice));

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
            at(*rep, to_size) = digit;
            to_size++;

            highest_power_rem /= base;
        }


        is_first = false;
    }


    //trim
    assert(to_size <= rep->size && "must be big enough");
    Char_Slice converted = {rep->data, to_size};

    size_t half_size = to_size / 2;
    for(size_t i = 0; i < half_size; i++)
    {
        char* l = &at(converted, i);
        char* r = &at(converted, to_size - i - 1);
        char temp = *l;
        *l = *r;
        *r = temp;
    }

    Char_Result result = {OK, converted};
    return result;
}

umax single_factorial_calculate(umax of_value)
{
    umax calculated = 1;
    for(umax factor = 2; factor <= of_value; factor++)
        calculated *= factor;

    return calculated;
}

umax single_factorial(umax of_value)
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

size_t required_factorial_out_size(umax of_value)
{
    if(of_value <= MAX_NATIVE_FACTORIAL)
        return 1;
    else
        return cast(size_t) (of_value - MAX_NATIVE_FACTORIAL);
}

Result factorial(Slice* out, umax of_value, const Optims* optims)
{
    if(out->size < required_factorial_out_size(of_value))
        return error_result(OUT_SIZE_TOO_SMALL);

    bool did_finish = true;
    if(of_value < MAX_NATIVE_FACTORIAL)
    {
        size_t fact = single_factorial(of_value);
        Slice output = from_number(out, fact);
        return ok_result(output);
    }

    const size_t highest_possible_native = single_factorial(MAX_NATIVE_FACTORIAL);
    Slice current_value = from_number(out, highest_possible_native);
    for(size_t factor = MAX_NATIVE_FACTORIAL + 1; factor <= of_value; factor++)
    {
        Result res = mul_short(out, to_cslice(current_value), cast(Digit) factor, optims);
        assert(res.state == OK && "shoudlnt fail");
        current_value = res.output;
    }

    return ok_result(current_value);
}