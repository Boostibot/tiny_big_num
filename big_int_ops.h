#pragma once

#include "preface.h"


#ifndef DO_RUNTIME_ONLY 
    //prevents compile time execution when on
    #define DO_RUNTIME_ONLY true
#endif

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
    #define MAX_RECURSION_DEPTH 10
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

struct Optim_Info
{
    bool mul_shift = DO_OPTIM_MUL_SHIFT;
    bool mul_consts = DO_OPTIM_MUL_CONSTS;
    bool mul_half_bits = DO_OPTIM_MUL_HALF_BITS;
    bool div_shift = DO_OPTIM_DIV_SHIFT;
    bool div_consts = DO_OPTIM_DIV_CONSTS;
    bool rem_optims = DO_OPTIM_REM_OPTIMS;

    size_t max_reusion_depth = MAX_RECURSION_DEPTH;
    size_t mul_quadratic_both_below_size = MUL_QUADRATIC_BOTH_BELOW_SIZE;
    size_t mul_quadratic_single_below_size = MUL_QUADRATIC_SINGLE_BELOW_SIZE;
    size_t trivial_pow_below_power = TRIVIAL_POW_BELOW_POWER;
};

template <typename T>
struct Trivial_Maybe
{
    bool has = false;
    T value = T();

    bool constexpr operator ==(Trivial_Maybe const&) const noexcept = default;
};

template <typename T>
func wrap(T const& value) -> Trivial_Maybe<T>
{
    return Trivial_Maybe<T>{true, value};
}

template <typename T>
func unwrap(Trivial_Maybe<T> const& maybe) -> T
{
    assert(maybe.has);
    return maybe.value;
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
func slice_size(Slice<T> slice, size_t from, size_t count)
{
    assert(from <= slice.size && from + count <= slice.size && "sliced portion must be entirely within the given slice");
    return Slice<T>{slice.data + from, count};
}

template <typename T>
func slice_range(Slice<T> slice, size_t from, size_t to)
{
    assert(to >= from);
    return ::slice_size<T>(slice, from, to - from);
}

template <typename T>
func slice(Slice<T> slice, size_t from) -> Slice<T> {
    return ::slice_size<T>(slice, from, slice.size - from);
}

template <typename T>
func trim(Slice<T> slice, size_t to_size) -> Slice<T> {
    return ::slice_size<T>(slice, 0, to_size);
}

template <typename T>
func are_aliasing(Slice<const T> left, Slice<const T> right) -> bool
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

    //const size_t diff = cast(size_t) abs(cast(ptrdiff_t) left.data - cast(ptrdiff_t) right.data);
    //return diff < left.size || diff < right.size;
}

template <typename T>
func are_one_way_aliasing(Slice<const T> before, Slice<const T> after) -> bool
{ 
    return (before.data + before.size > after.data) && (after.data > before.data);
}

template <typename T>
func check_are_aliasing(Slice<const T> left, Slice<const T> right) -> bool
{
    if constexpr(DO_RUNTIME_ONLY)
        return are_aliasing<T>(left, right);
    else
        return false;
}

template <typename T>
func check_are_one_way_aliasing(Slice<const T> before, Slice<const T> after) -> bool
{ 
    if constexpr(DO_RUNTIME_ONLY)
        return are_one_way_aliasing<T>(before, after);
    else
        return false;
}

enum class Iter_Direction 
{
    FORWARD,
    BACKWARD,
    NO_ALIAS
};
template <typename T>
proc copy_n(T* to, const T* from, size_t count, Iter_Direction direction) -> void
{
    Slice<const T> to_slice = {to, count};
    Slice<const T> from_slice = {from, count};

    if(direction == Iter_Direction::FORWARD)
        assert(check_are_one_way_aliasing(from_slice, to_slice) == false && "order must match aliasing");
    else if(direction == Iter_Direction::BACKWARD)
        assert(check_are_one_way_aliasing(to_slice, from_slice) == false && "order must match aliasing");
    else
        assert(check_are_aliasing(to_slice, from_slice) == false && "order must match aliasing");

    if(DO_RUNTIME_ONLY)
    {
        if(direction == Iter_Direction::NO_ALIAS)
            memcpy(to, from, count * sizeof(T));
        else
            memmove(to, from, count * sizeof(T));
        return;
    }

    if(direction == Iter_Direction::BACKWARD)
    {
        for(size_t i = count; i-- > 0;)
            to[i] = from[i];
    }
    else
    {
        for(size_t i = 0; i < count; i++)
            to[i] = from[i];
    }
}

template <typename T>
proc null_n(T* to, size_t count) -> void
{
    static_assert(std::is_integral_v<T>);
    if(DO_RUNTIME_ONLY)
    {
        memset(to, 0, count * sizeof(T));
        return;
    }

    for(size_t i = 0; i < count; i++)
        to[i] = 0;
}

template <typename T>
proc copy_slice(Slice<T>* to, Slice<const T> from, Iter_Direction direction) -> void
{
    assert(to->size == from.size && "sizes must match");
    return copy_n(to->data, from.data, from.size, direction);
}

template <typename T>
proc null_slice(Slice<T>* to)
{
    return null_n(to->data, to->size);
}


//@TODO: Hook to the existing parsing functions I have made
//@TODO: Make Better multiply algorhitm
//@TODO: Maybe just copy over 
//@TODO: Figure out param passing

template <typename T>
constexpr size_t BIT_SIZE = sizeof(T) * CHAR_BIT;

template <typename T>
constexpr size_t HALF_BIT_SIZE = BIT_SIZE<T> / 2;

template <typename T>
constexpr T FULL_MASK = cast(T)(-1);



template <typename T>
func set_bit(T field, size_t bit_pos, T val) -> T {
    assert(val == 0 || val == 1);
    assert(bit_pos < BIT_SIZE<T>);

    const T bit = val << bit_pos;
    const T mask = cast(T) 1 << bit_pos;
    const T res = (field & ~mask) | bit;
    return res;
};

template <typename T>
func get_bit(T field, size_t bit_pos) -> T 
{
    return (field >> bit_pos) & 1u;
};

template <typename T>
func set_nth_bit(Slice<T>* num, size_t i, T val)
{
    const size_t digit_i = i / BIT_SIZE<T>;
    const size_t bit_i = i % BIT_SIZE<T>;

    assert(digit_i < num->size);

    T* digit = &(*num)[digit_i];
    *digit = set_bit(*digit, bit_i, val);
};

template <typename T>
func get_nth_bit(Slice<T> num, size_t i) -> T {
    const size_t digit_i = i / BIT_SIZE<T>;
    const size_t bit_i = i % BIT_SIZE<T>;

    return get_bit(num[digit_i], bit_i);
};

template <typename T>
func high_mask(size_t index = HALF_BIT_SIZE<T>) -> T {
    assert(index < BIT_SIZE<T>);
    static_assert(std::is_unsigned_v<T>);
    return cast(T) (FULL_MASK<T> << index);
}

template <typename T>
func low_mask(size_t index = HALF_BIT_SIZE<T>) -> T {
    static_assert(std::is_unsigned_v<T>);
    return cast(T) ~high_mask<T>(index);
}

template <typename T>
func high_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
    static_assert(std::is_unsigned_v<T>);
    assert(index < BIT_SIZE<T>);
    return value >> index;
}

template <typename T>
func low_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
    static_assert(std::is_unsigned_v<T>);
    let mask = low_mask<T>(index);
    return value & mask;
};  

template <typename T>
func combine_bits(T low, T high, size_t index = HALF_BIT_SIZE<T>) -> T {
    assert(index < BIT_SIZE<T>);
    static_assert(std::is_unsigned_v<T>);
    return low_bits(low, index) | (high << index);
};

template <typename T>
func dirty_combine_bits(T low, T high, size_t index = HALF_BIT_SIZE<T>) -> T {
    assert(index < BIT_SIZE<T>);
    static_assert(std::is_unsigned_v<T>);
    assert(high_bits(low, index) == 0 && "low must not have high bits use combine_bits instead");
    return low | (high << index);
};

template <typename T>
func find_last_set_bit(T val) -> T  
{
    static_assert(std::is_unsigned_v<T>);
    if constexpr(sizeof(T) > 8)
    {
        //@TODO: upgrade to iterative on u64s and test
        T k = 0;

        while (val >>= 1) 
            k++;

        return k;
    }
    else
    {
        assert(sizeof(T) <= 8 && "Higher integer sizes not supported");
        T  k = 0;
        if (val > 0xFFFFFFFFu) { val >>= 32; k |= 32; }
        if (val > 0x0000FFFFu) { val >>= 16; k |= 16; }
        if (val > 0x000000FFu) { val >>= 8;  k |= 8;  }
        if (val > 0x0000000Fu) { val >>= 4;  k |= 4;  }
        if (val > 0x00000003u) { val >>= 2;  k |= 2;  }

        k |= (val & 2) >> 1;
        return k;
    }
}

template <typename T>
func find_first_set_bit(T val) -> T  
{
    if(val == 0)
        return 0;

    static_assert(std::is_unsigned_v<T>);
    if constexpr(sizeof(T) > 8)
    {
        T k = 0;

        while ((val & 1) == 0)
        {
            k++;
            val >>= 1;
        }
        return k;
    }
    else
    {
        //10110000
        //&   1111 -> nothing => +4 and shift
        //00001011
        //&     11 -> something => +0
        //00001011
        //&      1 -> somehing => +0
        // => answer: 4 + 0 + 0 = 4 

        assert(sizeof(T) <= 8 && "Higher integer sizes not supported");
        T  k = 0;

        if constexpr(sizeof(T) > 4)
            if ((val & 0xFFFFFFFFu) == 0) { val >>= 32; k |= 32; }
        if constexpr(sizeof(T) > 2)
            if ((val & 0x0000FFFFu) == 0) { val >>= 16; k |= 16; }
        if constexpr(sizeof(T) > 1)
            if ((val & 0x000000FFu) == 0) { val >>= 8;  k |= 8;  }

        if ((val & 0x0000000Fu) == 0) { val >>= 4;  k |= 4;  }
        if ((val & 0x00000003u) == 0) { val >>= 2;  k |= 2;  }

        k |= (~val & 0b1);
        return k;
    }
}

template <typename T>
func pop_count(T val) -> size_t
{
    static_assert(std::is_unsigned_v<T>);
    if constexpr(sizeof(T) < sizeof(uint32_t))
    {
        return pop_count<uint32_t>(cast(uint32_t) val);
    }
    else if constexpr(sizeof(T) == sizeof(uint32_t))
    {
        uint32_t i = cast(uint32_t) val;
        i = i - ((i >> 1) & 0x55555555);        // add pairs of bits
        i = (i & 0x33333333) + ((i >> 2) & 0x33333333);  // quads
        i = (i + (i >> 4)) & 0x0F0F0F0F;        // groups of 8
        return (i * 0x01010101) >> 24;          // horizontal sum of bytes
    }
    else if constexpr(sizeof(T) == sizeof(uint32_t)*2)
    {
        uint32_t low = cast(uint32_t) low_bits(val);
        uint32_t high = cast(uint32_t) high_bits(val);

        return pop_count<uint32_t>(low) + pop_count<uint32_t>(high);
    }
    else
    {
        //shoudl happen normally
        size_t count = 0;
        for (; val != 0; val >>= 1)
            if (val & 1)
                count++;

        return count;
    }
}


template <typename T>
struct Overflow
{
    T value;
    T overflow;

    bool constexpr operator ==(Overflow const&) const noexcept = default;
};

template <typename T>
struct Batch_Overflow
{
    Slice<T> slice;
    T overflow;

    bool constexpr operator ==(Batch_Overflow const&) const noexcept = default;
};

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

template <typename T>
func find_last_set_digit(Slice<T> num) -> size_t
{
    static_assert(std::is_unsigned_v<T>);

    size_t i = num.size;
    for(; i-- > 0;)
        if(num[i] != 0)
            break;

    return i;
}

template <typename T>
func find_first_set_digit(Slice<T> num) -> size_t
{
    static_assert(std::is_unsigned_v<T>);
    
    size_t i = 0;
    for(; i < num.size; i++)
        if(num[i] != 0)
            break;

    return i;
}

template <typename T>
func striped_trailing_zeros(Slice<T> num) -> Slice<T>
{
    return trim(num, find_last_set_digit(num) + 1);
}

template <typename T>
func is_striped_number(Slice<T> num) -> bool
{
    if(num.size == 0)
        return true;

    return last(num) != 0;
}

template <typename T>
func is_striped_representation(Slice<T> num) -> bool
{
    if(num.size == 0)
        return true;

    return num[0] != 0;
}

template <typename Digit, typename Number>
func digits_to_represent() -> size_t
{
    return (sizeof(Number) + sizeof(Digit) - 1) / sizeof(Digit);
}

template <typename To = std::uint64_t, typename T = size_t>
func to_number(Slice<T> bignum) -> Trivial_Maybe<To>
{
    static_assert(std::is_unsigned_v<T>);
    let stripped = striped_trailing_zeros(bignum);

    const size_t total_size = stripped.size * sizeof(T);
    if(total_size > sizeof(To))
        return {};

    To out = 0;
    for(size_t i = 0; i < stripped.size; i++)
    {
        size_t shift_by = i * BIT_SIZE<T>;
        assert(shift_by < BIT_SIZE<To>);
        out |= cast(To) stripped[i] << shift_by;
    }

    return wrap(out);
}

template <typename T, typename From = std::uint64_t>
proc from_number(Slice<T>* bignum, From from) -> Slice<T>
{
    static_assert(std::is_unsigned_v<T>);
    constexpr size_t count = digits_to_represent<T, From>();
    assert(bignum->size >= count);

    size_t i = 0;
    for(; i < count; i++)
    {
        From curr_digit = from >> (i * BIT_SIZE<T>);
        if(curr_digit == 0)
            break;

        (*bignum)[i] = cast(T) curr_digit;
    }

    return trim(*bignum, i);
}

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
//
//All this complexity could have been avoided if c++ had a proper standard instrinsics library.
template <typename T>
func add_overflow(T left, T right, T carry = 0) -> Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(high_bits(carry) == 0 && "carry must be half bits use add_overflow_any instead");
    const T res_low = low_bits(left) + low_bits(right) + carry;
    const T middle_carry = high_bits(res_low);

    const T res_high = high_bits(left) + high_bits(right) + middle_carry;
    const T new_carry = high_bits(res_high);

    return Overflow<T>{combine_bits(res_low, res_high), new_carry};
}

template <typename T>
func sub_overflow(T left, T right, T carry = 0) -> Overflow<T>
{
    assert(carry == 0 || carry == 1);

    const T res_low = low_bits(left) - low_bits(right) - carry;
    const T middle_carry = res_low >> (BIT_SIZE<T> - 1);

    const T res_high = high_bits(left) - high_bits(right) - middle_carry;
    const T new_carry = res_high >> (BIT_SIZE<T> - 1);

    return Overflow<T>{combine_bits(res_low, res_high), new_carry};
}

namespace detail
{
    template <typename T, typename... Ts>
    func add_overflow_any(Ts... args) -> Overflow<T>
    {
        static_assert((std::is_same_v<T, Ts> && ...));
        //will only overflow if we have more than 2^HALF_BIT_SIZE<T> + 1 operands => assert
        // (given that the smallest unsigned type is u8 we can safely sum up to 16 operands)
        // 
        //0x0F         <> 1 - from first
        //0x0F => 0x1E ^
        //0x0F => 0x2D |
        //0x0F => 0x3C | 0xF args == 2^HALF_BIT_SIZE<T> - 1
        //.....        |
        //0x0F => 0xF0 v
        //0x0F => 0xFF <> 1 
        //0x0F => 0x0E - !overflow!
        // => needs 2^HALF_BIT_SIZE<T> - 1 + 1 + 1 == 2^HALF_BIT_SIZE<T> + 1 operands
        constexpr size_t arg_count = sizeof...(args);
        constexpr size_t max_args = cast(size_t) (cast(T) 1 << HALF_BIT_SIZE<T>);
        static_assert(arg_count <= max_args && "too many operands");

        const T res_low = (low_bits(cast(T) args) + ...);
        const T middle_carry = high_bits(res_low);

        const T res_high = (high_bits(cast(T) args) + ...) + middle_carry;
        const T new_carry = high_bits(res_high);

        return Overflow<T>{combine_bits(res_low, res_high), new_carry};
    }

    template <typename T, typename... Ts>
    func sub_overflow_any(T left, Ts... args) -> Overflow<T>
    {
        static_assert((std::is_same_v<T, Ts> && ...));
        constexpr size_t arg_count = sizeof...(args);
        constexpr size_t max_args = cast(size_t) (cast(T) 1 << HALF_BIT_SIZE<T>);
        static_assert(arg_count < max_args && "too many operands");

        const T res_low = high_mask<T>() + low_bits(left) - (low_bits(args) + ...);
        const T middle_carry = low_mask<T>() - high_bits<T>(res_low);

        const T res_high = high_mask<T>() + high_bits(left) - (high_bits(args) + ...) - middle_carry;
        const T new_carry = low_mask<T>() - high_bits<T>(res_high);

        return Overflow<T>{combine_bits(res_low, res_high), new_carry};
    }

    template <typename T, typename... Ts>
    func add_no_overflow(Ts... args) -> T
    {
        static_assert((std::is_same_v<T, Ts> && ...));
        #ifdef NDEBUG
            return (args + ...);
        #else
            const size_t size = sizeof...(args);
            const T arr[] = {args...};
            T sum = 0;

            for(size_t i = 0; i < size; i++)
            {
                const T new_sum = sum + arr[i];
                assert(new_sum >= sum && "must not overflow");
                sum = new_sum;
            }

            return sum;
        #endif // NDEBUG
    }
}

template <typename T, typename... Ts>
func add_overflow_any(T left, T right, Ts... rest) -> Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    return detail::add_overflow_any<T>(left, right, rest...);
}

template <typename T, typename... Ts>
func sub_overflow_any(T left, T right, Ts... rest) -> Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    return detail::sub_overflow_any<T>(left, right, rest...);
}

template <typename T, typename... Ts>
func add_no_overflow(T left, T right, Ts... rest) -> T
{
    return detail::add_no_overflow<T>(left, right, rest...);
}

enum class Op_Location
{
    IN_PLACE,
    OUT_OF_PLACE
};

enum class Consume_Op
{
    ADD,
    SUB
};

template <typename T>
func consume_carry(Slice<T>* to, Slice<const T> left, T carry, size_t from, const Consume_Op consume_op, Op_Location location = Op_Location::OUT_OF_PLACE) -> Batch_Overflow<T>
{
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);

    const Slice<T> trimmed_to = trim(*to, left.size);
    size_t j = from;
    for (; carry != 0 && j < left.size != 0; j++)
    {
        const T digit = left[j];
        bool carry_consumed = false;
        T patch_res = 0;
        if(consume_op == Consume_Op::ADD)
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

    if(location == Op_Location::OUT_OF_PLACE)
        copy_n<T>(to->data + j, left.data + j, left.size - j, Iter_Direction::FORWARD);

    return Batch_Overflow<T>{trimmed_to, carry};
}


template <typename T>
func add_overflow_batch(Slice<T>* to, Slice<const T> left, T right, Op_Location location = Op_Location::OUT_OF_PLACE) -> Batch_Overflow<T>
{
    return consume_carry(to, left, right, 0, Consume_Op::ADD, location);
}

template <typename T>
func sub_overflow_batch(Slice<T>* to, Slice<const T> left, T right, Op_Location location = Op_Location::OUT_OF_PLACE) -> Batch_Overflow<T>
{
    return consume_carry(to, left, right, 0, Consume_Op::SUB, location);
}

template <typename T>
func add_overflow_batch(Slice<T>* to, Slice<const T> left, Slice<const T> right, Op_Location location = Op_Location::OUT_OF_PLACE, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= right.size);
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);
    assert(high_bits(carry_in) == 0);
    //assert(carry_in == 0 || carry_in == 1);

    if(left.size < right.size)
        std::swap(left, right);

    T carry = carry_in;
    for (size_t i = 0; i < right.size; i++)
    {
        let res = add_overflow<T>(left[i], right[i], carry);
        (*to)[i] = res.value;
        carry = res.overflow;
    }

    return consume_carry<T>(to, left, carry, right.size, Consume_Op::ADD, location);
}

template <typename T>
func sub_overflow_batch(Slice<T>* to, Slice<const T> left, Slice<const T> right, Op_Location location = Op_Location::OUT_OF_PLACE, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= right.size);
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);
    assert(carry_in == 0 || carry_in == 1);
    //assert(high_bits(carry_in) == 0);

    T carry = carry_in;
    proc update = [&](size_t i, Overflow<T> oveflow) {
        (*to)[i] = oveflow.value;
        carry = oveflow.overflow;
    };

    const size_t min_size = min(left.size, right.size);
    const size_t max_size = max(left.size, right.size);
    for (size_t i = 0; i < min_size; i++)
        update(i, sub_overflow<T>(left[i], right[i], carry));

    for (size_t i = left.size; i < right.size; i++)
        update(i, sub_overflow<T>(0, right[i], carry));

    if(right.size < left.size)
        return consume_carry<T>(to, left, carry, right.size, Consume_Op::SUB, location);

    return Batch_Overflow<T>{trim(*to, max_size), carry};
}

template <typename T>
func complement_overflow(T val, T carry) -> Overflow<T>
{
    const T val_inv = ~val;

    const T complement_low = low_bits(val_inv) + carry;
    const T middle_carry = high_bits(complement_low);

    const T complement_high = high_bits(val_inv) + middle_carry;
    const T new_carry = high_bits(complement_high);

    return Overflow<T>{combine_bits(complement_low, complement_high), new_carry};
}

template <typename T>
proc complement_overflow_batch(Slice<T>* to, Slice<const T> left, T carry_in = 1) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= left.size);
    assert(is_striped_number(left));
    assert(check_are_one_way_aliasing<T>(left, *to) == false);
    assert(carry_in == 1 || carry_in == 0);

    T carry = carry_in;
    for(size_t i = 0; i < left.size; i++)
    {
        let res = complement_overflow<T>(left[i], carry);
        (*to)[i] = res.value;
        carry = res.overflow;
    }

    return Batch_Overflow<T>{trim(*to, left.size), carry};
}

template <typename T>
func shift_up_overflow(T low_item, size_t by_bits, T high_item) -> Overflow<T>
{
    assert(by_bits < BIT_SIZE<T> && by_bits != 0 && "shift must be in valid range (cannot be zero because its unimplmentable without if)");
    const size_t remaining_bits = BIT_SIZE<T> - by_bits;

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

    const T low = high_bits<T>(low_item, remaining_bits);
    const T high = low_bits<T>(high_item, remaining_bits);

    const T composed = low | high << by_bits;

    const T out = high_bits(high_item, remaining_bits);
    const T shifted_out = out;
    return Overflow<T>{composed, shifted_out};
}

template <typename T>
func shift_down_overflow(T low_item, size_t by_bits, T high_item) -> Overflow<T>
{
    assert(by_bits < BIT_SIZE<T> && by_bits != 0);
    const size_t remaining_bits = BIT_SIZE<T> - by_bits;
    assert(remaining_bits < BIT_SIZE<T>);

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

    
    const T low = high_bits(low_item, by_bits);
    const T high = low_bits(high_item, by_bits);

    const T composed = low | high << remaining_bits;

    const T out = low_bits(low_item, by_bits);
    const T shifted_out = out << remaining_bits;
    return Overflow<T>{composed, shifted_out};
}

//We allow shifting while iterating in both directions:
// 
// FOR SHIFT_UP: the FORWARD direction is useful in general case
//  (for example when shifting array in place)
// and the BACKWARD direction can be used while shifting by more then one element
//  and writes to already passed memory are neccessary
// 
// FOR SHIFT_DOWN: its the polar opposite

//@TODO: make overflow consistent with division

template <typename T>
proc shift_up_overflow_batch(Slice<T>* out, Slice<const T> in, size_t by_bits, 
    Iter_Direction direction = Iter_Direction::FORWARD, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(out->size >= in.size);
    assert(by_bits < BIT_SIZE<T>);

    let out_trimmed = trim(*out, in.size);

    if(by_bits == 0)
    {
        copy_n<T>(out->data, in.data, in.size, direction);
        return Batch_Overflow<T>{out_trimmed, carry_in};
    }

    if(in.size == 0)
        return Batch_Overflow<T>{out_trimmed, carry_in};

    if(direction == Iter_Direction::FORWARD)
    {
        assert(check_are_one_way_aliasing<T>(in, *out) == false);
        T prev = carry_in;
        for (size_t i = 0; i < in.size; i++)
        {
            let curr = in[i];
            let res = shift_up_overflow<T>(prev, by_bits, curr);
            (*out)[i] = res.value;
            prev = curr;
        }

        let carry_out = shift_up_overflow<T>(prev, by_bits, 0);
        return Batch_Overflow<T>{out_trimmed, carry_out.value};
    }
    else
    {
        assert(check_are_one_way_aliasing<T>(*out, in) == false);
        const T shifted_out = in[in.size - 1];
        T prev = shifted_out;
        for (size_t i = in.size; i-- > 1; )
        {
            let curr = in[i - 1];
            (*out)[i] = shift_up_overflow<T>(curr, by_bits, prev).value;
            prev = curr;
        }

        {
            let res = shift_up_overflow<T>(carry_in, by_bits, prev);
            (*out)[0] = res.value;
        }


        let carry_out = shift_up_overflow<T>(shifted_out, by_bits, 0);
        return Batch_Overflow<T>{out_trimmed, carry_out.value};
    }
}

template <typename T>
proc shift_down_overflow_batch(Slice<T>* out, Slice<const T> in, size_t by_bits, 
    Iter_Direction direction = Iter_Direction::FORWARD, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(out->size >= in.size);
    assert(by_bits < BIT_SIZE<T>);

    let out_trimmed = trim(*out, in.size);
    if(by_bits == 0)
    {
        //@TODO: Handle carry_in! (consistently with the other edge case)
        copy_n<T>(out->data, in.data, in.size, direction);
        return Batch_Overflow<T>{out_trimmed, carry_in};
    }

    if(in.size == 0)
        return Batch_Overflow<T>{out_trimmed, carry_in};

    if(direction == Iter_Direction::FORWARD)
    {
        assert(check_are_one_way_aliasing<T>(in, *out) == false);

        const T shifted_out = in[0];
        T prev = shifted_out;
        for (size_t i = 0; i < in.size - 1; i++)
        {
            let curr = in[i + 1];
            let res = shift_down_overflow<T>(prev, by_bits, curr);
            (*out)[i] = res.value;
            prev = curr;
        }

        {
            let res = shift_down_overflow<T>(prev, by_bits, carry_in);
            (*out)[in.size - 1] = res.value;
        }

        let carry_out = shift_down_overflow<T>(0, by_bits, shifted_out);
        return Batch_Overflow<T>{out_trimmed, carry_out.value};
    }
    else
    {
        assert(check_are_one_way_aliasing<T>(*out, in) == false);

        T prev = carry_in;
        for (size_t i = in.size; i-- > 0;)
        {
            let curr = in[i];
            (*out)[i] = shift_down_overflow<T>(curr, by_bits, prev).value;
            prev = curr;
        }

        let carry_out = shift_down_overflow<T>(0, by_bits, prev);
        return Batch_Overflow<T>{out_trimmed, carry_out.value};
    }
}

//return -1 if left < right
//        0 if left == right
//        1 if left > right   
template <typename T>
func compare(Slice<const T> left, Slice<const T> right) -> int
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

template <typename T>
func is_equal(Slice<const T> left, Slice<const T> right) -> bool
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

template <typename T>
func mul_overflow(T left, T right, T last_value, Mul_Overflow_Optims const& mul_optims = Mul_Overflow_Optims::NONE) -> Overflow<T>
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
    
    size_t h = HALF_BIT_SIZE<T>;

    T l = left;
    T r = right;

    T l1 = low_bits(l); 
    T l2 = high_bits(l); 
    T r1 = low_bits(r);  
    T r2 = high_bits(r); 

    //calculate constants
    T s_cur = 0;
    T s_mixed_ns = 0;
    T s_mixed_ns_overflow = 0;
    T s_next_ns = 0;

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
        let s_mixed_ns_res = add_overflow<T>(l1*r2, l2*r1);

        s_mixed_ns = s_mixed_ns_res.value;
        s_mixed_ns_overflow = s_mixed_ns_res.overflow;
        s_cur = l1*r1;
        s_next_ns = l2*r2;
    }

    //split mixed
    assert(h < BIT_SIZE<T>);

    T low_mixed = low_bits(s_mixed_ns) << h; // even though we are taking lower half of s_mixed we have to insert it to upper half
    T high_mixed = high_bits(s_mixed_ns) + (s_mixed_ns_overflow << h); //add overflow as another digit

    //distribute mixed to appropriate halves

    let curr_value = add_overflow_any<T>(s_cur, low_mixed, last_value); //also add last_value to save ops
    T next_value =   add_no_overflow<T>(s_next_ns, high_mixed, curr_value.overflow); //also add the overflow

    return Overflow<T>{curr_value.value, next_value};
}
template <typename T>
func is_power_of_two(T num)
{
    return num != 0 && (num & (num - 1)) == 0;
}

template <typename T>
proc mul_overflow_batch(Slice<T>* to, Slice<const T> left, T right, Optim_Info const& optims, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);

    T carry = carry_in;
    let trimmed_to = trim(*to, left.size);

    if(optims.mul_consts)
    {
        if(left.size == 0)
            return Batch_Overflow<T>{trimmed_to, carry_in};

        if(right == 0) 
            return Batch_Overflow<T>{trim(*to, 0), 0};

        if(right == 1)
        {
            copy_n(to->data, left.data, left.size, Iter_Direction::FORWARD);
            return Batch_Overflow<T>{trimmed_to, 0};
        }
    }

    if(optims.mul_shift)
    {
        if(is_power_of_two(right))
        {
            //@NOTE: here carry is handled differently the in the rest of the mul algorhitm
            //  this doesnt affect how the algorhitm appears from the outside not does it affect the
            //  ability for batches to be chained together with carries but it forbids any kind of dependence
            //  on the carry in/out in this case. This is not ideal and should be changed so that the
            //  shift itself handles carry in similar way to mul
            let shift_bits = find_last_set_bit(right);
            return shift_up_overflow_batch<T>(to, left, shift_bits, Iter_Direction::FORWARD, carry_in);
        }
    }

    #define do_carry_loop(operation)            \
    {                                           \
        for (size_t i = 0; i < left.size; i++)  \
        {                                       \
            const Overflow<T> res = operation;  \
            (*to)[i] = res.value;               \
            carry = res.overflow;               \
        }                                       \
    }                                           \

    if(optims.mul_half_bits && high_bits(right) == 0)
        do_carry_loop((mul_overflow<T>(left[i], right, carry, Mul_Overflow_Optims::LOW_BITS_ONLY)))
    else if(optims.mul_half_bits && low_bits(right) == 0)
        do_carry_loop((mul_overflow<T>(left[i], right, carry, Mul_Overflow_Optims::HIGH_BITS_ONLY)))
    else
        do_carry_loop((mul_overflow<T>(left[i], right, carry)))

    #undef do_carry_loop
    return Batch_Overflow<T>{trimmed_to, carry};
}

template <typename T>
func div_overflow_low(T left, T right, T carry_in = 0) -> Overflow<T>
{
    //The algorhitm works as follows (only in different base - we use base 10 for demosntartion)
    // 61 / 5 == 10 + 11 / 5 == 10 + 2 == 12
    // 
    //To go through this instinctive algorhitm we did the followig:
    // {161/5} == [16/5]*10 + (16%5)*10 + {1/5} = 3*10 + {(10 + 1) / 5} == 30 + {11/5} == ... == 32 
    // where {} means normal (ie infinitely precise fraction) used for yet uncalculated division

    //this however only works for carry and right thats smaller than half bits 
    // (else it doesnt get carried through properly through the modulos
    //  so for example div_overflow_low(0, 0xFF, 0xF0) should return 0xF0 but returns less) => assert
    assert(high_bits(carry_in) == 0 && "carry must be single digit");
    assert(high_bits(right) == 0 && "right must be single digit");
    assert(right != 0 && "cannot divide by zero");

    const T operand_high = combine_bits<T>(high_bits(left), carry_in);
    const T res_high = operand_high / right;
    const T middle_carry = operand_high % right;

    const T operand_low = combine_bits<T>(low_bits(left), middle_carry);
    const T res_low = operand_low / right;
    const T out_carry = operand_low % right;

    const T res = dirty_combine_bits(res_low, res_high);
    return Overflow<T>{res, out_carry};
}


template <typename T>
proc div_overflow_low_batch(Slice<T>* to, Slice<const T> left, T right, Optim_Info const& optims, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= left.size);
    assert(high_bits(right) == 0 && "only works for divisors under half bit size");
    assert(right != 0 && "cannot divide by zero");
    assert(check_are_one_way_aliasing<T>(*to, left) == false);

    let trimmed_to = trim(*to, left.size);
    if(optims.div_consts)
    {
        if(right == 1)
        {
            copy_n(to->data, left.data, left.size, Iter_Direction::BACKWARD);
            return Batch_Overflow<T>{trimmed_to, 0}; 
        }
    }

    if(optims.div_shift)
    {
        if(is_power_of_two(right))
        {
            let shift_bits = find_last_set_bit(right);
            let shift_res = shift_down_overflow_batch<T>(to, left, shift_bits, Iter_Direction::BACKWARD, carry_in);
            //because of the way shifting down result is defined we have to process the reuslt
            const T carry = shift_res.overflow >> (BIT_SIZE<T> - shift_bits);
            return Batch_Overflow<T>{trimmed_to, carry};
        }
    }

    T carry = carry_in;
    for (size_t i = left.size; i-- > 0; )
    {
        let res = div_overflow_low<T>(left[i], right, carry);
        (*to)[i] = res.value;
        carry = res.overflow;
    }

    return Batch_Overflow<T>{trimmed_to, carry};
}

static constexpr bool DO_REM_OPTIMS = true;

template <typename T>
proc rem_overflow_low_batch(Slice<const T> left, T right, Optim_Info const& optims, T carry_in = 0) -> T
{
    static_assert(std::is_unsigned_v<T>);
    assert(high_bits(right) == 0 && "only works for divisors under half bit size");
    assert(right != 0 && "cannot divide by zero");

    if(optims.rem_optims)
    {
        if(is_power_of_two(right))
        {
            let shift_bits = find_last_set_bit(right);
            if(left.size == 0)
                return carry_in;

            return low_bits(left[0], shift_bits);
        }
    }

    T carry = carry_in;
    for (size_t i = left.size; i-- > 0; )
    {
        let res = div_overflow_low<T>(left[i], right, carry);
        carry = res.overflow;
    }

    return carry;
}


func required_add_to_size(size_t left_size, size_t right_size) -> size_t
{
    return max(left_size, right_size) + 1;
}

template <typename T>
proc add(Slice<T>* to, Slice<const T> left, Slice<const T> right, Op_Location location = Op_Location::OUT_OF_PLACE)
{
    static_assert(std::is_unsigned_v<T>);
    assert(is_striped_number(left));
    assert(is_striped_number(right));
    assert(to->size >= required_add_to_size(left.size, right.size));

    let res = add_overflow_batch<T>(to, left, right, location);
    if(res.overflow == 0)
    {
        if(location == Op_Location::IN_PLACE)
            return striped_trailing_zeros(res.slice);
        
        assert(is_striped_number(res.slice));
        return res.slice;
    }

    size_t slice_size = res.slice.size;
    Slice<T> trimmed_to = trim(*to, slice_size + 1);
    trimmed_to[slice_size] = res.overflow;

    assert(is_striped_number(trimmed_to));

    return trimmed_to;
}

func required_sub_to_size(size_t left_size, size_t right_size) -> size_t
{
    return left_size;
}

template <typename T>
proc sub(Slice<T>* to, Slice<const T> left, Slice<const T> right, Op_Location location = Op_Location::OUT_OF_PLACE)
{
    static_assert(std::is_unsigned_v<T>);
    assert(is_striped_number(left));
    assert(is_striped_number(right));
    assert(to->size >= required_sub_to_size(left.size, right.size));

    let res = sub_overflow_batch<T>(to, left, right, location);
    assert(res.overflow == 0 && "left must be bigger than right");

    return striped_trailing_zeros(res.slice);
}

func required_mul_to_size(size_t left_size, size_t right_size) -> size_t
{
    return left_size + right_size;
}

template <typename T>
proc mul(Slice<T>* to, Slice<const T> left, T right, Optim_Info const& optims)
{
    static_assert(std::is_unsigned_v<T>);
    assert(is_striped_number(left));
    assert(to->size >= required_sub_to_size(left.size, 1));

    let res = mul_overflow_batch<T>(to, left, right, optims);
    if(res.overflow == 0)
    {
        assert(is_striped_number(res.slice));
        return res.slice;
    }

    size_t slice_size = res.slice.size;
    Slice<T> trimmed_to = trim(*to, slice_size + 1);
    trimmed_to[slice_size] = res.overflow;

    assert(is_striped_number(trimmed_to));
    return trimmed_to;
}

template <typename T>
struct Div_Res
{
    Slice<T> quotient;
    Slice<T> remainder;
};


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

template <typename T, bool DO_QUOTIENT = true>
proc div_bit_by_bit(Slice<T>* quotient, Slice<T>* remainder, Slice<const T> num, Slice<const T> den, Optim_Info const& optims) -> Trivial_Maybe<Div_Res<T>>
{
    static_assert(std::is_unsigned_v<T>);
    assert(is_striped_number(num));
    assert(is_striped_number(den));
    assert(check_are_aliasing<T>(*quotient, num) == false);
    assert(check_are_aliasing<T>(*remainder, num) == false);
    assert(check_are_aliasing<T>(*quotient, den) == false);
    assert(check_are_aliasing<T>(*remainder, den) == false);

    const size_t required_quot_size = DO_QUOTIENT ? required_div_quotient_size(num.size, den.size) : 0;
    const size_t required_rem_size = required_div_remainder_size(num.size, den.size);

    if(DO_QUOTIENT)
        assert(check_are_aliasing<T>(*quotient, *remainder) == false);

    if(den.size == 0)
        return {};

    if(num.size == 0)
    {
        let stripped_quotient = trim(*quotient, 0);
        let stripped_remainder = trim(*remainder, 0);
        return wrap(Div_Res<T>{stripped_quotient, stripped_remainder});
    }

    if(num.size < den.size)
    {
        assert(remainder->size >= num.size);
        mut stripped_remainder = trim(*remainder, num.size);
        mut stripped_quotient = trim(*quotient, 0);

        copy_slice<T>(&stripped_remainder, num, Iter_Direction::NO_ALIAS);
        return wrap(Div_Res<T>{stripped_quotient, stripped_remainder});
    }

    //@TODO: move up
    Slice<T> trimmed_remainder = trim(*remainder, required_rem_size);
    Slice<T> trimmed_quotient = trim(*quotient, required_quot_size);    

    if(DO_QUOTIENT)
        null_slice(&trimmed_quotient);

    null_slice(&trimmed_remainder);
    if(den.size == 1 && high_bits(last(den)) == 0)
    {
        assert(trimmed_remainder.size > 0 && "at this point should not be 0");

        Slice<T> stripped_quotient;
        T remainder_val;
        if(DO_QUOTIENT)
        {
            let div_res = div_overflow_low_batch<T>(&trimmed_quotient, num, last(den), optims);
            stripped_quotient = striped_trailing_zeros<T>(div_res.slice);
            remainder_val = div_res.overflow;
        }
        else
        {
            remainder_val = rem_overflow_low_batch<T>(num, last(den), optims);
            stripped_quotient = trim(trimmed_quotient, 0);
        }

        size_t remainder_size = 0;
        if(remainder_val != 0)
        {
            trimmed_remainder[0] = remainder_val;
            remainder_size = 1;
        }

        let stripped_remainder = trim(trimmed_remainder, remainder_size);
        return wrap(Div_Res<T>{stripped_quotient, stripped_remainder});
    }

    Slice<T> curr_remainder = trim(trimmed_remainder, 0);
    for(size_t i = BIT_SIZE<T> * num.size; i-- > 0; ) 
    {
        let last_size = curr_remainder.size;
        let shift_res = shift_up_overflow_batch<T>(&curr_remainder, curr_remainder, 1);
        const T num_ith_bit = get_nth_bit(num, i);

        if(shift_res.overflow != 0 || (last_size == 0 && num_ith_bit == 1))
        {
            curr_remainder = trim(trimmed_remainder, last_size + 1);
            curr_remainder[last_size] = shift_res.overflow;
        }

        set_nth_bit(&trimmed_remainder, 0, num_ith_bit);
        if(compare<T>(curr_remainder, den) >= 0)
        {
            let sub_res = sub_overflow_batch<T>(&curr_remainder, curr_remainder, den, Op_Location::IN_PLACE);
            curr_remainder = striped_trailing_zeros<T>(curr_remainder);
            assert(sub_res.overflow == 0 && "should not overflow");
                
            if(DO_QUOTIENT)
                set_nth_bit<T>(&trimmed_quotient, i, 1);
        }
    }

    assert(is_striped_number(curr_remainder));
    let stripped_quotient = DO_QUOTIENT 
        ? striped_trailing_zeros<T>(trimmed_quotient)
        : trim(trimmed_quotient, 0);

    return wrap(Div_Res<T>{stripped_quotient, curr_remainder});
}

template <typename T>
proc div(Slice<T>* quotient, Slice<T>* remainder, Slice<const T> num, Slice<const T> den, Optim_Info const& optims) -> Trivial_Maybe<Div_Res<T>>
{
    return div_bit_by_bit(quotient, remainder, num, den, optims);
}

template <typename T>
proc rem_bit_by_bit(Slice<T>* remainder, Slice<const T> num, Slice<const T> den, Optim_Info const& optims) -> Trivial_Maybe<Slice<T>>
{
    static_assert(std::is_unsigned_v<T>);
    Slice<T> quotient = {nullptr, 0};
    let div_res = div_bit_by_bit<T, false>(&quotient, remainder, num, den, optims);
    if(div_res.has == false)
        return {};

    assert(div_res.quotient.size == 0 && "in case of only remainder the quotient size should be 0");
    assert(div_res.quotient.data == nullptr && "and data nullptr as we set it");
    return wrap(div_res.remainder);
}

template <typename T>
proc fused_mul_add_overflow_batch(Slice<T>* to, Slice<const T> added, Slice<const T> multiplied, T coeficient, 
    Optim_Info const& optims, const Op_Location location = Op_Location::OUT_OF_PLACE, T add_carry = 0, T mul_carry = 0) -> Batch_Overflow<T>
{   
    assert(added.size >= multiplied.size);
    assert(to->size >= added.size);
    assert(check_are_one_way_aliasing<T>(added, *to) == false);
    assert(check_are_one_way_aliasing<T>(multiplied, *to) == false);

    if(optims.mul_consts)
    {
        const Slice<T> trimmed_to = trim(*to, added.size);
        if(coeficient == 0) 
        {
            if (location == Op_Location::IN_PLACE)
                return Batch_Overflow<T>{trim(*to, 0), 0};
            else
            {
                copy_n(to->data, added.data, added.size, Iter_Direction::FORWARD);
                return Batch_Overflow<T>{trimmed_to, 0};
            }
        }

        if(coeficient == 1)
            return add_overflow_batch<T>(to, added, multiplied, location);
    }

    Slice<T> trimmed_to = trim(*to, added.size);
    if(optims.mul_shift 
        && optims.mul_consts //zero check must be performed (else shifting is UB!)
        && is_power_of_two(coeficient))
    {
        let shift_bits = find_last_set_bit(coeficient);
        const Slice<T> trimmed_to = trim(*to, added.size);
        T prev_muled = mul_carry;
        for (size_t j = 0; j < multiplied.size; j++)
        {
            const T curr_muled = multiplied[j];
            const T curr_added = added[j];

            let mul_res = shift_up_overflow<T>(prev_muled, shift_bits, curr_muled);
            let add_res = add_overflow<T>(curr_added, mul_res.value, add_carry);

            trimmed_to[j] = add_res.value;

            prev_muled = curr_muled;
            add_carry = add_res.overflow;
        }

        mul_carry = shift_up_overflow<T>(prev_muled, shift_bits, 0).value;
    }
    else
    {
        for (size_t j = 0; j < multiplied.size; j++)
        {
            T curr_muled = multiplied[j];
            T curr_added = added[j];

            let mul_res = mul_overflow<T>(curr_muled, coeficient, mul_carry);
            let add_res = add_overflow<T>(curr_added, mul_res.value, add_carry);

            trimmed_to[j] = add_res.value;

            mul_carry = mul_res.overflow;
            add_carry = add_res.overflow;
        }
    }
    
    const T combined_carry = add_no_overflow(mul_carry, add_carry);
    return consume_carry<T>(&trimmed_to, added, combined_carry, multiplied.size, Consume_Op::ADD, location);
}

template <typename T>
proc mul(Slice<T>* to, Slice<T>* aux, Slice<const T> left, Slice<const T> right, Optim_Info const& optims, size_t depth = 0);

func required_mul_quadratic_auxiliary_size(size_t left_size, size_t right_size) -> size_t
{
    return required_mul_to_size(left_size, right_size);
}

template <typename T>
proc mul_quadratic(Slice<T>* to, Slice<T>* temp, Slice<const T> left, Slice<const T> right, Optim_Info const& optims) -> Slice<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(check_are_aliasing<T>(*to, left) == false);
    assert(check_are_aliasing<T>(*to, right) == false);
    assert(check_are_aliasing<T>(*temp, left) == false);
    assert(check_are_aliasing<T>(*temp, right) == false);
    assert(check_are_aliasing<T>(*to, *temp) == false);

    if(left.size < right.size)
        std::swap(left, right);

    const size_t max_result_size = required_mul_to_size(left.size, right.size);
    const size_t max_temp_size = required_mul_to_size(left.size, right.size);
    assert(to->size >= max_result_size);
    assert(temp->size >= max_temp_size);
    
    Slice<T> trimmed_to = trim(*to, max_result_size);
    Slice<T> trimmed_temp = trim(*temp, max_temp_size);
    null_n(trimmed_to.data, trimmed_to.size);

    for(size_t i = 0; i < right.size; i++)
    {
        let mul_res = mul_overflow_batch<T>(&trimmed_temp, left, right[i], optims);

        let temp_num = unwrap(to_number<u64>(mul_res.slice));

        //mul_overflow_batch is limited to left.size elems => in case of overflow
        // even if the overflow would fit into temp it is not added => we add it back
        let mul_size = mul_res.slice.size;
        let mul_slice = trim(trimmed_temp, mul_size + 1); 
        mul_slice[mul_size] = mul_res.overflow;

        mut shifted_to = slice(trimmed_to, i);
        assert(shifted_to.size >= mul_slice.size);

        let add_res = add_overflow_batch<T>(&shifted_to, shifted_to, mul_slice, Op_Location::IN_PLACE);
        assert(add_res.overflow == 0 && "no final carry should be left");
    }

    return striped_trailing_zeros(trimmed_to);
}

//@TODO: make it so that to is only operated on to the required sized
template <typename T>
proc mul_quadratic_fused(Slice<T>* to, Slice<const T> left, Slice<const T> right, Optim_Info const& optims) -> Slice<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(check_are_aliasing<T>(*to, left) == false);
    assert(check_are_aliasing<T>(*to, right) == false);

    if(left.size < right.size)
        std::swap(left, right);

    const size_t max_result_size = required_mul_to_size(left.size, right.size);
    assert(to->size >= max_result_size);
    Slice<T> trimmed_to = trim(*to, max_result_size);
    null_n(trimmed_to.data, trimmed_to.size);

    //Slice<T> last_modified = slice(trimmed_to, 0); //@TODO
    for(size_t i = 0; i < right.size; i++)
    {
        mut shifted_to = slice(trimmed_to, i);

        let mul_add_res = fused_mul_add_overflow_batch<T>(&shifted_to, shifted_to, left, right[i], optims, Op_Location::IN_PLACE);

        assert(mul_add_res.overflow == 0 && "no carry should happen in this case");
    }

    return striped_trailing_zeros(trimmed_to);
}


func required_mul_karatsuba_auxiliary_size(size_t left_size, size_t right_size) -> size_t
{
    if(left_size < right_size)
        std::swap(left_size, right_size);

    size_t base = max(left_size, right_size) / 2; 
    size_t capped_digits = min(right_size, base);

    size_t x1 = base;
    size_t x2 = left_size - base;

    size_t y1 = capped_digits;
    size_t y2 = right_size - capped_digits;

    size_t required_add_x_sum_size = required_add_to_size(x1, x2);
    size_t required_add_y_sum_size = required_add_to_size(y1, y2);

    size_t required_z2_size = required_mul_to_size(required_add_x_sum_size, required_add_y_sum_size);

    return required_z2_size 
        + required_add_x_sum_size 
        + required_add_y_sum_size;
}


template <typename T>
proc mul_karatsuba(Slice<T>* to, Slice<T>* aux, Slice<const T> left, Slice<const T> right, Optim_Info const& optims, size_t depth = 0, bool is_run_alone = true) -> Slice<T>
{
    static_assert(std::is_unsigned_v<T>);

    assert(is_striped_number(left));
    assert(is_striped_number(right));

    //@TODO make function for this
    assert(are_aliasing<T>(*to, left) == false);
    assert(are_aliasing<T>(*aux, left) == false);
    assert(are_aliasing<T>(*to, right) == false);
    assert(are_aliasing<T>(*aux, right) == false);
    assert(are_aliasing<T>(*aux, *to) == false);

    Slice<const T> x = left;
    Slice<const T> y = right;

    //Perfomrs the karatsuba multiplicaton step
    //x*y = (x1*b + x2)(y1*b + y2) = (x1*y1)b^2 + (x1*y2 + x2*y1)b + (x2*y2)
    //    = z1*b^2 + z2*b + z3

    // z1 = x1*y1
    // z3 = x2*y2
    // z2 = x1*y2 + x2*y1 = (x1 + x2)(y1 + y2) - z1 - z3

    if(x.size < y.size)
        std::swap(x, y);

    if(is_run_alone)
    {
        if(x.size == 0 || y.size == 0)
            return trim(*to, 0);

        if(y.size == 1)
            return mul<T>(to, x, y[0], optims);
    }

    //size_t split_digits = min(x.size / 2, y.size / 2);
    //we split according to the bigger one. This might seem odd but in the case where we will split too much
    // and y1 will be 0 we will have one less multiplication to worry about
    size_t base = max(x.size, y.size) / 2; 
    size_t capped_digits = min(y.size, base);
    assert(base != 0);

    Slice<const T> x1 = slice(x, base);
    Slice<const T> x2 = striped_trailing_zeros(trim(x, base));

    Slice<const T> y1 = slice(y, capped_digits);
    Slice<const T> y2 = striped_trailing_zeros(trim(y, capped_digits));


    size_t required_out_size = required_mul_to_size(x.size, y.size);
    Slice<T> trimmed_out = trim(*to, required_out_size);

    size_t required_add_x_sum_size = required_add_to_size(x1.size, x2.size);
    size_t required_add_y_sum_size = required_add_to_size(y1.size, y2.size);

    size_t required_z1_size = required_mul_to_size(x1.size, y1.size);
    size_t required_z2_size = required_mul_to_size(required_add_x_sum_size, required_add_y_sum_size); //after multiplies will be even less
    size_t required_z3_size = required_mul_to_size(x2.size, y2.size);

    size_t total_used_aux_size = required_add_x_sum_size + required_add_y_sum_size + required_z2_size;
    size_t returned_required_aux = required_mul_karatsuba_auxiliary_size(left.size, right.size);
    assert(total_used_aux_size <= returned_required_aux);

    Slice<T> remaining_aux = slice(*aux, total_used_aux_size);

    size_t z1_from_index = trimmed_out.size - required_z1_size;

    //we place z1 and z3 directly into the output buffer 
    // into appropiate positions as if multiplied by base
    Slice<T> z1_slot = slice(trimmed_out, z1_from_index);
    Slice<T> z3_slot = trim(trimmed_out, required_z3_size);

    assert(are_aliasing<T>(z1_slot, z3_slot) == false);

    //we multyply into z3 and null the area between it and begginign of z1 
    // (ie the place where z2 will get added to)
    Slice<T> z3 = mul<T>(&z3_slot, &remaining_aux, x2, y2, optims, depth + 1);
    Slice<T> z3_to_z1 = slice_range(trimmed_out, z3.size, z1_from_index);
    null_slice(&z3_to_z1);

    //we multiply into z1 and null the remaining size to the end of trimmed_out
    // this way trimmed out should look like the following:
    // 000[ z1 ]000[ z3 ]
    Slice<T> z1 = mul<T>(&z1_slot, &remaining_aux, x1, y1, optims, depth + 1);
    Slice<T> z1_up = slice(trimmed_out, z1.size + z1_from_index);
    null_slice(&z1_up);

    size_t used_to = 0;
    Slice<T> z2_slot = slice_size(*aux, used_to, required_z2_size);
    used_to += required_z2_size;

    Slice<T> x_sum_slot = slice_size(*aux, used_to, required_add_x_sum_size);
    used_to += required_add_x_sum_size;

    Slice<T> y_sum_slot = slice_size(*aux, used_to, required_add_y_sum_size);
    used_to += required_add_y_sum_size;

    Slice<T> x_sum = add<T>(&x_sum_slot, x1, x2);
    Slice<T> y_sum = add<T>(&y_sum_slot, y1, y2);

    Slice<T> x_sum_y_sum = mul<T>(&z2_slot, &remaining_aux, x_sum, y_sum, optims, depth + 1);
    Slice<T> x_sum_y_sum_m_z1 = sub<T>(&x_sum_y_sum, x_sum_y_sum, z1, Op_Location::IN_PLACE);
    Slice<T> z2 = sub<T>(&x_sum_y_sum_m_z1, x_sum_y_sum_m_z1, z3, Op_Location::IN_PLACE);

    //instead of multiplying z2 by base we add it to the appropriate position
    Slice<T> out_z2_up = slice(trimmed_out, base);

    let add_res = add_overflow_batch<T>(&out_z2_up, out_z2_up, z2, Op_Location::IN_PLACE);
    assert(add_res.overflow == 0 && "should not overflow");

    Slice<T> out = striped_trailing_zeros(trimmed_out);

    return out;
}

template <typename T>
proc mul(Slice<T>* to, Slice<T>* aux, Slice<const T> left, Slice<const T> right, Optim_Info const& optims, size_t depth)
{
    static_assert(std::is_unsigned_v<T>);
    assert(is_striped_number(left));
    assert(is_striped_number(right));
    assert(are_aliasing<T>(*to, left) == false);
    assert(are_aliasing<T>(*to, right) == false);
    assert(to->size >= required_mul_to_size(left.size, right.size));

    Slice<T> ret;
   
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
        && depth < optims.max_reusion_depth)
    {
        size_t required_aux = required_mul_karatsuba_auxiliary_size(max_size, min_size);
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

template <typename T>
func log2(Slice<const T> left) -> size_t
{
    assert(is_striped_number(left));
    if(left.size == 0)
        return 0;

    size_t last_log = find_last_set_bit(last(left));
    return last_log + (left.size - 1) * BIT_SIZE<T>;
}

func required_pow_to_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
{
    if(power == 0 || num_bit_size == 0)
        return 1;

    size_t bit_size = (num_bit_size + 1) * power;
    size_t item_size = div_round_up(bit_size, single_digit_bit_size) + 1;

    return item_size;
}

func required_pow_by_squaring_single_auxiliary_swap_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
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

func required_pow_by_squaring_auxiliary_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
{
    if(power == 0)
        return 0;
    
    size_t square = required_pow_by_squaring_single_auxiliary_swap_size(num_bit_size, power, single_digit_bit_size);
    size_t out_swap = required_pow_to_size(num_bit_size, power, single_digit_bit_size);

    return square * 2 + out_swap;
}


template <typename T>
func required_pow_to_size(Slice<const T> num, size_t power) -> size_t
{
    return required_pow_to_size(log2(num), power, BIT_SIZE<T>);
}
template <typename T>
func required_pow_by_squaring_single_auxiliary_swap_size(Slice<const T> num, size_t power) -> size_t
{
    return required_pow_by_squaring_single_auxiliary_swap_size(log2(num), power, BIT_SIZE<T>);
}
template <typename T>
func required_pow_by_squaring_auxiliary_size(Slice<const T> num, size_t power) -> size_t
{
    return required_pow_by_squaring_auxiliary_size(log2(num), power, BIT_SIZE<T>);
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

template <typename T>
proc pow_by_squaring(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t power, Optim_Info const& optims) -> Slice<T>
{
    assert(is_striped_number(num));
    //no two can alias
    assert(are_aliasing<T>(*to, num) == false);
    assert(are_aliasing<T>(*aux, num) == false);
    assert(are_aliasing<T>(*aux, *to) == false);

    const size_t required_size = required_pow_to_size(num, power);
    const size_t required_sigle_aux = required_pow_by_squaring_single_auxiliary_swap_size<T>(num, power);

    assert(to->size >= required_size);
    assert(aux->size >= required_pow_by_squaring_auxiliary_size(num, power));

    if(power == 0 || (num.size == 1 && num[0] == 1))
    {
        mut out = trim(*to, 1);
        out[0] = 1;
        return out;
    }

    if(num.size == 0)
        return trim(*to, 0);

    if(optims.mul_consts)
    {
        if(power == 1)
        {
            copy_n(to->data, num.data, num.size, Iter_Direction::NO_ALIAS);
            return trim(*to, num.size);
        }
        if(power == 2)
            return mul<T>(to, aux, num, num, optims);
    }

    Slice<T> output_aux1 = *to;
    Slice<T> output_aux2 = trim(*aux, required_size);
    
    Slice<T> after_out_aux = slice(*aux, required_size);

    Slice<T> square_aux1 = slice_size(after_out_aux, required_sigle_aux*0, required_sigle_aux);
    Slice<T> square_aux2 = slice_size(after_out_aux, required_sigle_aux*1, required_sigle_aux);

    Slice<T> remianing_aux = slice(after_out_aux, 2*required_sigle_aux);

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
        std::swap(output_aux1, output_aux2);

    Slice<const T> curr_square = num;
    Slice<T> curr_output = trim(output_aux1, 1);
    curr_output[0] = 1;

    for(size_t i = 0; i <= max_pos; i++)
    {
        size_t bit = get_bit(power, i);
        if(bit)
        {
            //curr_output *= curr_square;
            curr_output = mul<T>(&output_aux2, &remianing_aux, curr_output, curr_square, optims);
            std::swap(output_aux1, output_aux2);
        }

        //if we are gonna stop next iteration dont waste time
        if(i != max_pos)
        {
            //curr_square *= curr_square;
            curr_square = mul<T>(&square_aux2, &remianing_aux, curr_square, curr_square, optims);
            std::swap(square_aux1, square_aux2);
        }
    }

    assert(curr_output.data == to->data);
    assert(is_striped_number(curr_output));
    return curr_output;
}

template <typename T>
func required_trivial_pow_auxiliary_size(Slice<const T> num, size_t power) -> size_t
{
    if(power <= 1)
        return 0;

    //minus one since in the algorhitm we always multiply from output buffer to auxiliary
    // so that the final mutliply will be to the output buffer
    // => the maximum power that will be stored in the auxiliary is one less
    return required_pow_to_size(num, power - 1);
}


template <typename T>
proc trivial_pow(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t power, Optim_Info const& optims) -> Slice<T>
{
    assert(is_striped_number(num));

    assert(are_aliasing<T>(*to, num) == false);
    assert(are_aliasing<T>(*aux, num) == false);
    assert(are_aliasing<T>(*aux, *to) == false);

    const size_t required_to_size = required_pow_to_size(num, power);
    const size_t required_aux_size = required_trivial_pow_auxiliary_size(num, power);
    assert(to->size >= required_to_size);
    assert(aux->size >= required_aux_size);

    if(power == 0 || (num.size == 1 && num[0] == 1))
    {
        mut out = trim(*to, 1);
        out[0] = 1;
        return out;
    }

    if(num.size == 0)
        return trim(*to, 0);

    Slice<T> output_aux1 = trim(*to, required_to_size);
    Slice<T> output_aux2 = trim(*aux, required_aux_size);

    Slice<T> remianing_aux = slice(*aux, required_aux_size);

    // we shoudl start with an appropiate one so that we finnish in the *to buffer
    const size_t num_assignments = power - 1;
    if(num_assignments % 2 != 0)
        std::swap(output_aux1, output_aux2);

    Slice<T> curr = trim(output_aux1, num.size);
    copy_slice<T>(&curr, num, Iter_Direction::NO_ALIAS);

    for(size_t i = 0; i < num_assignments; i++)
    {
        Slice<T> next = mul<T>(&output_aux2, &remianing_aux, curr, num, optims);
        std::swap(output_aux1, output_aux2);
        curr = next;
    }

    assert(curr.data == to->data && "we should finish in the output buffer");
    assert(is_striped_number(curr));
    return curr;
}

template <typename T>
proc pow(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t power, Optim_Info const& optims) -> Slice<T>
{
    size_t at_least_to = required_pow_to_size(num, power);
    size_t at_least_aux = required_trivial_pow_auxiliary_size(num, power);
    assert(to->size >= at_least_to);
    assert(aux->size >= at_least_aux);

    if(power < optims.trivial_pow_below_power)
        return trivial_pow<T>(to, aux, num, power, optims);

    const size_t required_to_size = required_pow_to_size(num, power);
    const size_t required_aux_size = required_pow_by_squaring_auxiliary_size(num, power);

    if(to->size >= required_to_size && aux->size >= required_aux_size)
        return pow_by_squaring(to, aux, num, power, optims);

    return trivial_pow<T>(to, aux, num, power, optims);

    const size_t bits = (log2(num) + 1);

    //upper estimate for single swap aux
    //size_t required_sing_swap_aux_size_rough = (bits*power / BIT_SIZE<T> + 2);
    //size_t total_aux_size_rough = (bits*power / BIT_SIZE<T> + 2) * 2;
    //=> 
    // required = ([bits*power / BIT_SIZE<T>] + 2) * 2
    // required/2 - 2 = [bits*power / BIT_SIZE<T>]
    // (required/2 - 2) * BIT_SIZE<T> = [bits*power] - bits*power%BIT_SIZE<T>
    // [bits*power] = (required/2 - 2) * BIT_SIZE<T> + bits*power%BIT_SIZE<T>
    // power = ((required/2 - 2) * BIT_SIZE<T> + bits*power%BIT_SIZE<T>) / bits
    // power < (required/2 - 2) * BIT_SIZE<T> / bits

    size_t max_aux_power = (aux->size/2 - 2) * BIT_SIZE<T> / bits;
    size_t max_to_power = (aux->size - 2) * BIT_SIZE<T> / bits;
   
    size_t max_pow_by_squaring_power = min(max_aux_power, max_to_power);
    size_t remianing_power = power - max_pow_by_squaring_power;
    //@TODO: try one power higher than the one caluclated and see if that fit
    // (rounding errors)

    Slice<T> half_powed = pow_by_squaring(to, aux, num, max_pow_by_squaring_power, optims);
    Slice<T> new_num = trim(*aux, half_powed.size);
    Slice<T> remaining_aux = slice(*aux, half_powed.size);
    copy_slice<T>(&new_num, *to, Iter_Direction::NO_ALIAS);

    return trivial_pow<T>(to, &remaining_aux, new_num, remianing_power, optims);
}


func root_estimate_bit_size(size_t num_bit_size, size_t root) -> size_t {
    //const u64 upper_initial_estimate = cast(u64) 1 << ((log2x + n) / n);
    size_t estimate_bit_size = (num_bit_size + root) / root;
    return estimate_bit_size;
}

func root_required_to_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size) -> size_t {
    if(root == 0)
        return 1;

    if(root == 1)
        return div_round_up(num_bit_size, single_digit_bit_size);

    //the +1 is for the addition of the quotient we perform in place
    return div_round_up(root_estimate_bit_size(num_bit_size, root), single_digit_bit_size) + 1;
}

func root_required_aux_size(size_t num_bit_size, size_t root, size_t single_digit_bit_size) -> size_t {

    size_t other_estimate_size = root_required_to_size(num_bit_size, root, single_digit_bit_size);
    size_t pow_aux_size = required_pow_by_squaring_auxiliary_size(num_bit_size, root, single_digit_bit_size);
    size_t pow_to_size = required_pow_by_squaring_auxiliary_size(num_bit_size, root, single_digit_bit_size);

    return other_estimate_size + pow_aux_size + pow_to_size;
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

    const u64 upper_initial_estimate = cast(u64) 1 << ((log2x + n) / n);
    u64 r = upper_initial_estimate;
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


template <typename T>
proc root(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t root, Optim_Info const& optims) -> Slice<T>
{
    size_t num_bit_size = log2(num);   
    size_t required_to_size = root_required_to_size(num.size, root, BIT_SIZE<T>);
    size_t required_aux_size = root_required_aux_size(num_bit_size, root, BIT_SIZE<T>);

    assert(to->size >= required_to_size);
    assert(aux->size >= required_aux_size);

    assert(high_bits(cast(T) root) == 0 && "only small roots allowed - would require generalization of this algorhitm");
    if(num.size == 0)
    {
        to->data[0] = 1;
        return trim(*to, 1);
    }

    if(root <= 0)
    {
        Slice<T> ret = trim(*to, num.size);
        copy_slice<T>(&ret, num, Iter_Direction::NO_ALIAS);
        return ret;
    }

    if(num.size == 1)
    {
        to->data[0] = cast(T) iroot_shifting(cast(u64) num[0], cast(u64) root);
        return trim(*to, 1);
    }
    size_t estimate_bit_size = root_estimate_bit_size(num_bit_size, root);

    Slice<T> out1 = *to;
    Slice<T> out2 = trim(*aux, required_to_size);
    Slice<T> rest_aux = slice(*aux, required_to_size);

    Slice<T> initial_estimate = out2;

    null_slice(&initial_estimate);
    set_nth_bit<T>(&initial_estimate, estimate_bit_size, 1);
    initial_estimate = striped_trailing_zeros(initial_estimate);
    //assert(is_striped_number(curr_estimate));

    Slice<T> curr_estimate = initial_estimate;
    Slice<T> prev_estimate = {};


    while(true)
    {
        prev_estimate = curr_estimate;
        assert(is_striped_number(prev_estimate));

        //new_upper_r = (root-1) * prev_estimate + num / ipow_squares(prev_estimate, root-1);

        
        //the sizes we will need to store the results:
        size_t powed_size = required_pow_to_size<T>(prev_estimate, root - 1);
        size_t muled_size = required_mul_to_size(prev_estimate.size, 1);

        //assert(powed_size <= num.size && "is it?");
        
        //we will perform the operations in order of most memory consumption.
        // that way they can use as much still-not-occupied storage to speed up their execution
        Slice<T> pow_to = trim(rest_aux, powed_size);
        Slice<T> after_pow_aux = slice(rest_aux, powed_size);
        Slice<T> powed = pow<T>(&pow_to, &after_pow_aux, prev_estimate, root - 1, optims);

        size_t div_size = required_div_quotient_size(num.size, powed.size) + 1;
        size_t rem_size = required_div_remainder_size(num.size, powed.size);

        Slice<T> div_to = trim(after_pow_aux, div_size);
        Slice<T> after_div_aux = slice(after_pow_aux, div_size);
        Slice<T> rem_to = trim(after_div_aux, rem_size);


        Div_Res<T> res = unwrap(div<T>(&div_to, &rem_to, num, powed, optims));
        Slice<T> dived = res.quotient;
        assert(is_striped_number(dived));

        size_t add_size = required_add_to_size(muled_size, dived.size);

        Slice<T> new_estimate_to = out1;
        assert(new_estimate_to.data != prev_estimate.data && "storages must differ");
        assert(new_estimate_to.size >= add_size || true &&  
            "should fit into the prepared slot - maybe will fail because add_size is theoretical reuquired size"
            "but new_estimate_to.size is the actual real size which might - and probably be - smaller than the theoretical size"
            "we might have to lie about size to the algorhitm here!");

        //we perform the multiply add fused
        Batch_Overflow<T> fused_result;
        
        if(dived.size != 0)
            fused_result = fused_mul_add_overflow_batch<T>(&new_estimate_to, dived, prev_estimate, cast(T) root - 1, optims);
        else
            fused_result = mul_overflow_batch<T>(&new_estimate_to, prev_estimate, cast(T) root - 1, optims);
            
        assert(fused_result.overflow == 0 && "should not overflow");
        assert(is_striped_number(fused_result.slice));
        std::swap(out1, out2);


        //r = new_upper_r / n;
        Slice<T> new_upper_estimate = fused_result.slice;
        Slice<T> new_estimate_pre = div_overflow_low_batch<T>(&new_upper_estimate, new_upper_estimate, cast(T) root, optims).slice;
        Slice<T> new_estimate = striped_trailing_zeros(new_estimate_pre);
        curr_estimate = new_estimate;

        //if(r >= prev_r)
            //break; //and the result is in prev_estimate (estimates shoudl be equal but just to be safe)
        int compared = compare<T>(new_estimate, prev_estimate);
        if(compared >= 0)
            break;
    }

    if(prev_estimate.data == to->data)
        return prev_estimate;

    //if the current output is in the wrong storage copy over
    Slice<T> dest = trim(*to, prev_estimate.size);
    copy_slice<T>(&dest, prev_estimate, Iter_Direction::NO_ALIAS);
    return dest;
}


template <typename T>
proc reverse(Slice<T>* arr) -> void
{
    Slice<T> ref = *arr;
    size_t half_size = ref.size / 2;
    for(size_t i = 0; i < half_size; i++)
        std::swap(ref[i], ref[ref.size - i - 1]);
}

runtime_func size_to_fit_base(double from_base, double from_size, double to_base) -> double
{
    //const double max_from = pow(from_base, from_size);
    //const double to_size = log(max_from) / log(to_base);
    assert(from_base > 1 && to_base > 1 && "base must be bigger than 2");
    const double to_size = from_size * log(from_base) / log(to_base);
    return to_size;
}

//what is the size of the output buffer when converting from number to its in base rep?
template <typename Num>
runtime_func required_size_to_base(size_t num_size, size_t to_base) -> size_t
{
    double in_one_digit = round(pow(2, BIT_SIZE<Num>));
    return cast(size_t) ceil(size_to_fit_base(in_one_digit, cast(double) num_size, cast(double) to_base));
}

//what is the size of the output buffer when converting from base rep to number?
template <typename Num>
runtime_func required_size_from_base(size_t represenation_size, size_t from_base) -> size_t
{
    //to_size = from_size * log(from_base) / log(to_base);
    double bit_size = cast(double) represenation_size * log(cast(double) from_base) / log(2.0);
    size_t max_bits = cast(size_t) ceil(bit_size);
    size_t digits = (max_bits + BIT_SIZE<Num> - 1) / BIT_SIZE<Num>;
    return digits;
}


struct To_Base_State
{
    size_t index = 0;
    size_t buffer_size = 0;
};

template <typename Rep>
struct To_Base_Result
{
    Rep value;
    To_Base_State state;
    bool finished = true;
};


template <typename Num>
func to_base_init(Slice<const Num> num) -> To_Base_State
{
    assert(is_striped_number(num));
    return To_Base_State{0, num.size};
}

template <typename Num>
proc to_base_convert(To_Base_State state, Slice<Num>* temp, Slice<const Num> num, Num base, Optim_Info const& optims) -> To_Base_Result<Num>
{
    static_assert(std::is_unsigned_v<Num>);
    assert(is_striped_number(num));
    assert(temp->size >= num.size && "temp must be big enough");
    assert(base >= 2 && "base must be bigger than two");
    assert(high_bits(base) == 0 && "base must be low bits only of Num type (so that short div algorhimt can be used)");

    if(state.buffer_size == 0)
        return To_Base_Result<Num>{0, state, true};

    size_t buff_size = state.buffer_size;
    Slice<Num> curr_buff = trim(*temp, buff_size);
    Slice<const Num> curr_div_from;
    if(state.index == 0)  
        curr_div_from = trim(num, buff_size);
    else
        curr_div_from = curr_buff;

    let res = div_overflow_low_batch<Num>(&curr_buff, curr_div_from, base, optims);
    let digit = res.overflow;
    if(last(curr_buff) == 0)
    {
        assert(buff_size != 0);
        buff_size = buff_size - 1;
    }

    let new_state = To_Base_State{state.index + 1, buff_size};
    return To_Base_Result<Num>{digit, new_state, false};
}

template <typename Rep>
proc to_base_finish(Slice<Rep>* rep) -> Slice<Rep>
{
    static_assert(std::is_unsigned_v<Rep>);
    reverse(rep);
    return *rep;
}

template <typename Num, typename Rep, typename Conversion_Fn>
proc to_base(Slice<Rep>* rep, Slice<Num>* temp, Slice<const Num> num, Num base, Conversion_Fn conversion, Optim_Info const& optims) -> Slice<Rep>
{
    static_assert(std::is_unsigned_v<Num>);
    static_assert(sizeof(Rep) <= sizeof(Num));

    using Converted = decltype(conversion(cast(Num) 0));
    static_assert(std::is_same_v<Converted, Rep>);


    if constexpr(std::is_same_v<Num, Rep>)
        assert(check_are_aliasing<Num>(num, *rep) == false);

    assert(is_striped_number(num));
    assert(temp->size >= num.size && "temp must be big enough");
    assert(base >= 2 && "base must be bigger than two");
    assert(high_bits(base) == 0 && "base must be low bits only of Num type (so that short div algorhimt can be used)");

    if constexpr (DO_RUNTIME_ONLY)
    {
        size_t at_least = required_size_to_base<Num>(num.size, base);
        assert(at_least <= rep->size && "there must be enough size rep represent the number even in the worst case scenario");
    }

    bool is_first = true;
    size_t to_size = 0;

    if(num.size == 0)
        return trim(*rep, 0);

    //convert
    Slice<Num> val_left = *temp;
    while(true)
    {
        const Slice<const Num> curr_div_from = is_first ? num : val_left;

        if(curr_div_from.size == 0)
            break;

        let res = div_overflow_low_batch<Num>(&val_left, curr_div_from, base, optims);
        val_left = striped_trailing_zeros(val_left);
        let digit = conversion(res.overflow);
        (*rep)[to_size] = digit;
        to_size++;

        is_first = false;
    }


    //finish
    Slice<Rep> converted = trim(*rep, to_size);
    assert(is_striped_number(converted));
    reverse(&converted);
    return converted;
}


struct From_Base_State
{
    size_t index = 0;
    size_t buffer_size = 0;
    size_t required_size = 0;
};

func from_base_init() -> From_Base_State
{
    return {0, 0, 1};
}

template <typename Num>
proc from_base_convert(From_Base_State state, Slice<Num>* num, Num digit, Num base, Optim_Info const& optims) -> From_Base_State
{
    static_assert(std::is_unsigned_v<Num>);

    //digit could have also had in general bigger type than Num
    // (we would just convert from the number to static buffer and then 
    //  add the result) but I believe that is such a specific case it would
    //  just make testing harder 
    assert(base > digit && "all digits must be valid");

    size_t curr_size = state.buffer_size;
    size_t required_size = curr_size + 1;

    if(num->size < required_size)
    {
        From_Base_State new_state = state;
        new_state.required_size = required_size;
        return new_state;
    }
    
    Slice<Num> curr_num = trim(*num, curr_size);
    size_t overflows_had = 0;

    let mul_res = mul_overflow_batch<Num>(&curr_num, curr_num, base, optims);
    if(mul_res.overflow != 0)
    {
        curr_size += 1;
        overflows_had += 1;
        curr_num = trim(*num, curr_size);
        curr_num[curr_size - 1] = mul_res.overflow;
    }

    let add_res = add_overflow_batch<Num>(&curr_num, curr_num, digit, Op_Location::IN_PLACE);
    if(add_res.overflow != 0)
    {
        curr_size += 1;
        overflows_had += 1;
        curr_num = trim(*num, curr_size);
        curr_num[curr_size - 1] = add_res.overflow;
    }

    assert(overflows_had < 2 && "it shouldnt be possible to overflow twice when the digit is under base");
    From_Base_State new_state = state;
    new_state.index += 1;
    new_state.buffer_size = curr_size;

    return new_state;
}

template <typename Num>
proc from_base_finish(Slice<Num>* num) -> Slice<Num>
{
    static_assert(std::is_unsigned_v<Num>);
    return striped_trailing_zeros(*num);
}

/*
template <typename Num, typename Rep>
proc from_base2(Slice<Num>* num, Slice<const Rep> rep, Num base, Optim_Info const& optims) -> Slice<Num>
{
    static_assert(std::is_unsigned_v<Num>);
    static_assert(sizeof(Num) >= sizeof(Rep));
    assert(is_striped_representation(rep));
    if constexpr(std::is_same_v<Rep, Num>)
        assert(check_are_aliasing<Rep>(rep, *num) == false);

    if constexpr (DO_RUNTIME_ONLY)
    {
        size_t at_least = required_size_from_base<Num>(rep.size, base);
        assert(at_least <= num->size && "there must be enough size in num to represent the number even in the worst case scenario");
    }

    null_n(num->data, num->size);
    mut state = from_base_init();

    while(state.index < rep.size)
    {
        assert(state.required_size <= num->size);
        Num digit = cast(Num) rep[state.index];
        let new_state = from_base_convert<Num>(state, num, digit, base, optims);
        state = new_state;
    }

    Slice<Num> trimmed_num = trim(*num, state.buffer_size);
    return from_base_finish(&trimmed_num);
}

template <typename Num, typename Rep>
proc to_base2(Slice<Rep>* rep, Slice<Num>* temp, Slice<const Num> num, Num base, Optim_Info const& optims) -> Slice<Rep>
{
    static_assert(std::is_unsigned_v<Num>);
    static_assert(sizeof(Rep) >= sizeof(Num));
    if constexpr(std::is_same_v<Num, Rep>)
        assert(check_are_aliasing<Num>(num, *rep) == false);

    if constexpr (DO_RUNTIME_ONLY)
    {
        size_t at_least = required_size_to_base<Num>(num.size, base);
        assert(at_least <= rep->size && "there must be enough size rep represent the number even in the worst case scenario");
    }

    To_Base_State state = to_base_init(num);
    while(true)
    {
        let res = to_base_convert(state, temp, num, base, optims);
        if(res.finished)
            break;

        (*rep)[state.index] = cast(Rep) res.value;
        state = res.state;
    }

    Slice<Rep> trimmed_rep = trim(*rep, state.index);
    return to_base_finish(&trimmed_rep);
}
*/

//100_000_000

template <typename Num, typename Rep, typename Conversion_Fn>
proc from_base(Slice<Num>* num, Slice<const Rep> rep, Num base, Conversion_Fn conversion, Optim_Info const& optims) -> Slice<Num>
{
    static_assert(std::is_unsigned_v<Num>);
    static_assert(sizeof(Num) >= sizeof(Rep));

    using Converted = decltype(conversion(cast(Rep) 0));
    static_assert(std::is_same_v<Converted, Num>);

    assert(is_striped_representation(rep));
    if constexpr(std::is_same_v<Rep, Num>)
        assert(check_are_aliasing<Rep>(rep, *num) == false);

    if constexpr (DO_RUNTIME_ONLY)
    {
        size_t at_least = required_size_from_base<Num>(rep.size, base);
        assert(at_least <= num->size && "there must be enough size num represent the number even in the worst case scenario");
    }

    null_n(num->data, num->size);

    size_t from_size = sizeof(Rep);
    size_t to_size = sizeof(Num);

    for(size_t i = 0; i < rep.size; i++)
    {
        Num digit = conversion(rep[i]);
        assert(base > digit && "all digits must be valid");

        Slice<Num> added = {&digit, 1};
        let mul_res = mul_overflow_batch<Num>(num, *num, base, optims);
        let add_res = add_overflow_batch<Num>(num, *num, digit, Op_Location::IN_PLACE);
        //let add_res = add_overflow_batch<Num>(num, *num, added, Op_Location::IN_PLACE);

        assert(add_res.overflow == 0);
        assert(mul_res.overflow == 0);
    }

    return striped_trailing_zeros(*num);
}