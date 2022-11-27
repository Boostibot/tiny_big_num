#pragma once

#include <bit>
#include "preface.h"

using Max_Unsigned_Type = u64;

static constexpr bool DO_DEFAULT_MUL_SHIFT_OPTIMS = false;
static constexpr bool DO_DEFAULT_MUL_CONST_OPTIMS = false;
static constexpr bool DO_DEFAULT_MUL_HALF_BITS_OPTIMS = false;

static constexpr bool DO_DEFAULT_DIV_SHIFT_OPTIMS = false;
static constexpr bool DO_DEFAULT_DIV_CONST_OPTIMS = false;
static constexpr bool DO_DEFAULT_DIV_FUSED_SHIFT_OPTIMS = true;
static constexpr bool DO_RUNTIME_ONLY_CHECKS = true; //prevents compile time execution when on

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

    //bool constexpr operator ==(Slice const&) const noexcept = default;
    constexpr operator Slice<const T>() const noexcept
    {
        return Slice<const T>{this->data, this->size};
    }
};

template <typename T>
func unsafe_slice(Slice<T> slice, size_t from, size_t count) -> Slice<T>
{
    return Slice<T>{slice.data + from, count};
}

template <typename T>
func unsafe_slice(Slice<T> slice, size_t from) -> Slice<T>
{
    return Slice<T>{slice.data + from, slice.size - from};
}

template <typename T>
func slice(Slice<T> slice, size_t from, size_t count)
{
    assert(from <= slice.size && from + count <= slice.size && "sliced portion must be entirely within the given slice");
    return unsafe_slice<T>(slice, from, count);
}

template <typename T>
func slice(Slice<T> slice, size_t from) -> Slice<T>
{
    return ::slice<T>(slice, from, slice.size - from);
}

template <typename T>
func trim(Slice<T> slice, size_t to_size) -> Slice<T>
{
    return ::slice<T>(slice, 0, to_size);
}

template <typename T>
func are_aliasing(Slice<const T> left, Slice<const T> right) -> bool
{ 
    const size_t diff = cast(size_t) abs(cast(ptrdiff_t) left.data - cast(ptrdiff_t) right.data);
    return diff < left.size || diff < right.size;
}

template <typename T>
func are_one_way_aliasing(Slice<const T> before, Slice<const T> after) -> bool
{ 
    return (before.data + before.size > after.data) && (after.data > before.data);
}

template <typename T>
func check_are_aliasing(Slice<const T> left, Slice<const T> right) -> bool
{
    if constexpr(DO_RUNTIME_ONLY_CHECKS)
        return are_aliasing<T>(left, right);
    else
        return false;
}

template <typename T>
func check_are_one_way_aliasing(Slice<const T> before, Slice<const T> after) -> bool
{ 
    if constexpr(DO_RUNTIME_ONLY_CHECKS)
        return are_one_way_aliasing<T>(before, after);
    else
        return false;
}

enum class Iter_Direction 
{
    FORWARD,
    BACKWARD
};
template <typename T>
proc copy_n(T* to, const T* from, size_t count, Iter_Direction direction) -> void
{
    Slice<const T> to_slice = {to, count};
    Slice<const T> from_slice = {from, count};
    if(direction == Iter_Direction::FORWARD)
    {
        assert(check_are_one_way_aliasing(from_slice, to_slice) == false && "order must match aliasing");
        for(size_t i = 0; i < count; i++)
            to[i] = from[i];
    }
    else
    {
        assert(check_are_one_way_aliasing(to_slice, from_slice) == false && "order must match aliasing");
        for(size_t i = count; i-- > 0;)
            to[i] = from[i];
    }
}

template <typename T>
proc null_n(T* to, size_t count) -> void
{
    static_assert(std::is_integral_v<T>);
    for(size_t i = 0; i < count; i++)
        to[i] = 0;
}

//@TODO: Make optim struct and default optim struct that will get passed during runtime
//        (for individual ops require compile time individual switches)
//        (this will most likely when just const passed get optimized out)
//@TODO: Figure out how to check if multiplied or divided number is just low batch mutiplication but shifted (how to count the first bit?)
//@TODO: Add tests that test slightly (off by one) aliasing instead of full overlap to determine if asserts work correctly
//@TODO: Tidy up the div proc - move bit manipulation out
//@TODO: Merge typed and utyped tests so that the typed part executes only when the supplied typed matches
//@TODO: Figure out the smallest possible requirements for div
//@TODO: Insatntiate different optims
//@TODO: Create interface for parsing and printing functions (taking functions probably or returnign state to be executed repeatedly) 
//@TODO: Hook to the existing parsing functions I have made
//@TODO: Make Better multiply algorhitm
//@TODO: Maybe just copy over 
//@TODO: Add define flags setting default
//@TODO: Figure out param passing

template <typename T>
constexpr size_t BIT_SIZE = sizeof(T) * CHAR_BIT;

template <typename T>
constexpr size_t HALF_BIT_SIZE = BIT_SIZE<T> / 2;

template <typename T>
constexpr T FULL_MASK = cast(T) Max_Unsigned_Type(-1);


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
func is_striped_form(Slice<T> num) -> bool
{
    if(num.size == 0)
        return true;

    return last(num) != 0;
}

template <typename T>
func is_reverse_striped_form(Slice<T> num) -> bool
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

template <typename To = Max_Unsigned_Type, typename T = size_t>
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

template <typename T, typename From = Max_Unsigned_Type>
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
    //we separate it like so that the user must supply at least 2 operands
    return detail::add_overflow_any<T>(left, right, rest...);
}

template <typename T, typename... Ts>
func add_no_overflow(T left, T right, Ts... rest) -> T
{
    return detail::add_no_overflow<T>(left, right, rest...);
}

template <typename T>
func sub_overflow(T left, T right, T carry = 0) -> Overflow<T>
{
    const T res_low = low_bits(left) - low_bits(right) - carry;
    const T middle_carry = res_low >> (BIT_SIZE<T> - 1);

    const T res_high = high_bits(left) - high_bits(right) - middle_carry;
    const T new_carry = res_high >> (BIT_SIZE<T> - 1);

    return Overflow<T>{combine_bits(res_low, res_high), new_carry};
}

template <typename T, bool IN_PLACE = false>
func add_overflow_batch_patch(Slice<T>* to, Slice<const T> left, T carry, size_t from) -> Batch_Overflow<T>
{
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);

    const Slice<T> trimmed_to = trim(*to, left.size);
    size_t j = from;
    if(carry != 0)
    {
        for (; j < left.size != 0; j++)
        {
            const T digit = left[j];
            const T patch_res = digit + carry;
            trimmed_to[j] = patch_res;

            if(patch_res > digit)
            {
                carry = 0;
                j++;
                break;
            }
        }
    }

    if constexpr (IN_PLACE == false)
        copy_n<T>(to->data + j, left.data + j, left.size - j, Iter_Direction::FORWARD);

    return Batch_Overflow<T>{trimmed_to, carry};
}

template <typename T, bool IN_PLACE = false>
func sub_overflow_batch_patch(Slice<T>* to, Slice<const T> left, T carry, size_t from) -> Batch_Overflow<T>
{
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);

    const Slice<T> trimmed_to = trim(*to, left.size);
    size_t j = from;
    if(carry != 0)
    {
        for (; j < left.size != 0; j++)
        {
            const T digit = left[j];
            const T patch_res = digit - carry;
            trimmed_to[j] = patch_res;

            if(patch_res < digit)
            {
                carry = 0;
                j++;
                break;
            }
        }
    }

    if constexpr (IN_PLACE == false)
        copy_n<T>(to->data + j, left.data + j, left.size - j, Iter_Direction::FORWARD);

    return Batch_Overflow<T>{trimmed_to, carry};
}

template <typename T, bool IN_PLACE = false>
func add_overflow_batch(Slice<T>* to, Slice<const T> left, Slice<const T> right, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= right.size);
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);
    assert(high_bits(carry_in) == 0);

    if(left.size < right.size)
        std::swap(left, right);

    T carry = carry_in;
    for (size_t i = 0; i < right.size; i++)
    {
        let res = add_overflow<T>(left[i], right[i], carry);
        (*to)[i] = res.value;
        carry = res.overflow;
    }

    return add_overflow_batch_patch<T, IN_PLACE>(to, left, carry, right.size);
}

template <typename T, bool IN_PLACE = false>
func sub_overflow_batch(Slice<T>* to, Slice<const T> left, Slice<const T> right, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= right.size);
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);
    assert(high_bits(carry_in) == 0);

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
        return sub_overflow_batch_patch<T, IN_PLACE>(to, left, carry, right.size);

    //for (size_t i = right.size; i < left.size; i++)
        //update(i, sub_carry_patch_step<T>(left[i], carry));

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
    assert(is_striped_form(left));
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

    const T out = low_bits(high_item, remaining_bits);
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
template <typename T>
proc shift_up_overflow_batch(Slice<T>* out, Slice<const T> in, size_t by_bits, 
    T carry_in = 0, Iter_Direction direction = Iter_Direction::FORWARD) -> Batch_Overflow<T>
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
proc shift_down_overflow_batch(Slice<T>* out, Slice<const T> in, size_t by_bits, T carry_in = 0, 
    Iter_Direction direction = Iter_Direction::FORWARD) -> Batch_Overflow<T>
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
    assert(is_striped_form(left));
    assert(is_striped_form(right));
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

template <typename T, Mul_Overflow_Optims OPTIMS = Mul_Overflow_Optims::NONE>
func mul_overflow(T left, T right, T last_value = 0) -> Overflow<T>
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
    if constexpr(OPTIMS == Mul_Overflow_Optims::LOW_BITS_ONLY)
    {
        s_cur = l1*r1;
        s_mixed_ns = l2*r1;
    }
    else if(OPTIMS == Mul_Overflow_Optims::HIGH_BITS_ONLY)
    {
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

template <typename T, 
    bool DO_SHIFT = DO_DEFAULT_MUL_SHIFT_OPTIMS, 
    bool DO_CONSTS = DO_DEFAULT_MUL_CONST_OPTIMS, 
    bool DO_HALVES = DO_DEFAULT_MUL_HALF_BITS_OPTIMS>
proc mul_overflow_batch(Slice<T>* to, Slice<const T> left, T right, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= left.size);
    assert(check_are_one_way_aliasing<T>(left, *to) == false);

    T carry = carry_in;
    let trimmed_to = trim(*to, left.size);

    //@TODO: compress options - make detail function out of loop
    if constexpr(DO_CONSTS)
    {
        if(right == 0) 
            return Batch_Overflow<T>{trim(*to, 0), 0};

        if(right == 1)
        {
            copy_n(to->data, left.data, left.size, Iter_Direction::FORWARD);
            return Batch_Overflow<T>{trimmed_to, 0};
        }
    }

    if constexpr(DO_SHIFT)
    {
        if(is_power_of_two(right))
        {
            let shift_bits = find_last_set_bit(right);
            return shift_up_overflow_batch<T>(to, left, shift_bits, carry_in, Iter_Direction::FORWARD);
        }
    }

    if constexpr(DO_HALVES)
    {
        if(high_bits(right) == 0)
        {
            for (size_t i = 0; i < left.size; i++)
            {
                let res = mul_overflow<T, Mul_Overflow_Optims::LOW_BITS_ONLY>(left[i], right, carry);
                (*to)[i] = res.value;
                carry = res.overflow;
            }
            return Batch_Overflow<T>{trimmed_to, carry};
        }

        if(low_bits(right) == 0)
        {
            for (size_t i = 0; i < left.size; i++)
            {
                let res = mul_overflow<T, Mul_Overflow_Optims::HIGH_BITS_ONLY>(left[i], right, carry);
                (*to)[i] = res.value;
                carry = res.overflow;
            }
            return Batch_Overflow<T>{trimmed_to, carry};
        }
    }

    for (size_t i = 0; i < left.size; i++)
    {
        let res = mul_overflow<T>(left[i], right, carry);
        (*to)[i] = res.value;
        carry = res.overflow;
    }

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


template <typename T, 
    bool DO_SHIFT = DO_DEFAULT_DIV_SHIFT_OPTIMS, 
    bool DO_CONST = DO_DEFAULT_DIV_CONST_OPTIMS>
proc div_overflow_low_batch(Slice<T>* to, Slice<const T> left, T right, T carry_in = 0) -> Batch_Overflow<T>
{
    static_assert(std::is_unsigned_v<T>);
    assert(to->size >= left.size);
    assert(high_bits(right) == 0 && "only works for divisors under half bit size");
    assert(right != 0 && "cannot divide by zero");
    assert(check_are_one_way_aliasing<T>(*to, left) == false);

    let trimmed_to = trim(*to, left.size);
    if constexpr(DO_CONST)
    {
        if(right == 1)
        {
            copy_n(to->data, left.data, left.size, Iter_Direction::BACKWARD);
            return Batch_Overflow<T>{trimmed_to, 0}; 
        }
    }

    if constexpr(DO_SHIFT)
    {
        if(is_power_of_two(right))
        {
            let shift_bits = find_last_set_bit(right);
            let shift_res = shift_down_overflow_batch<T>(to, left, shift_bits, carry_in, Iter_Direction::BACKWARD);
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

template <typename T, bool DO_OPTIMS = DO_REM_OPTIMS>
proc rem_overflow_low_batch(Slice<const T> left, T right, T carry_in = 0) -> T
{
    static_assert(std::is_unsigned_v<T>);
    assert(high_bits(right) == 0 && "only works for divisors under half bit size");
    assert(right != 0 && "cannot divide by zero");

    if constexpr(DO_OPTIMS)
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


template <typename T>
struct Div_Res
{
    Slice<T> quotient;
    Slice<T> remainder;
};

template <typename T, 
    bool DO_SHIFT = DO_DEFAULT_DIV_SHIFT_OPTIMS, 
    bool DO_CONST = DO_DEFAULT_DIV_CONST_OPTIMS, 
    bool DO_QUOTIENT = true>
proc div_bit_by_bit(Slice<T>* quotient, Slice<T>* remainder, Slice<const T> num, Slice<const T> den) -> Trivial_Maybe<Div_Res<T>>
{
    static_assert(std::is_unsigned_v<T>);
    assert(is_striped_form(num));
    assert(is_striped_form(den));
    assert(check_are_aliasing<T>(*quotient, num) == false);
    assert(check_are_aliasing<T>(*remainder, num) == false);
    assert(check_are_aliasing<T>(*quotient, den) == false);
    assert(check_are_aliasing<T>(*remainder, den) == false);

    if(DO_QUOTIENT)
        assert(check_are_aliasing<T>(*quotient, *remainder) == false);

    if(den.size == 0)
        return {};

    const size_t min_size = min(num.size, den.size);
    assert(remainder->size >= min_size);
    if(DO_QUOTIENT)
    {
        assert(quotient->size >= min_size);
        null_n(quotient->data, quotient->size);
    }

    null_n(remainder->data, remainder->size);

    if(num.size < den.size)
    {
        copy_n(remainder->data, num.data, num.size, Iter_Direction::FORWARD);

        let stripped_quotient = trim(*quotient, 0);
        let stripped_remainder = trim(*remainder, num.size);
        return wrap(Div_Res<T>{stripped_quotient, stripped_remainder});
    }

    if(den.size == 1 && high_bits(last(den)) == 0)
    {
        if constexpr (DO_QUOTIENT == false)
        {
            const T rem_res = rem_overflow_low_batch<T, DO_SHIFT, DO_CONST>(num, last(den));
            remainder[0] = rem_res;

            let stripped_quotient = trim(*quotient, 0);
            let stripped_remainder = trim(*remainder, 1);
            return wrap(Div_Res<T>{stripped_quotient, stripped_remainder});
        }
        else
        {
            let div_res = div_overflow_low_batch<T, DO_SHIFT, DO_CONST>(quotient, num, last(den));
            assert(remainder->size > 0 && "at this point should not be 0");

            size_t remainder_size = 0;
            if(div_res.overflow != 0)
            {
                (*remainder)[0] = div_res.overflow;
                remainder_size = 1;
            }

            let stripped_quotient = striped_trailing_zeros<T>(div_res.slice);
            let stripped_remainder = trim(*remainder, remainder_size);
            return wrap(Div_Res<T>{stripped_quotient, stripped_remainder});
        }
    }

    constexpr let set_bit = [](T field, size_t bit_pos, T val) -> T {
        assert(val == 0 || val == 1);
        assert(bit_pos < BIT_SIZE<T>);

        const T bit = val << bit_pos;
        const T mask = cast(T) 1 << bit_pos;
        const T res = (field & ~mask) | bit;
        return res;
    };

    constexpr let get_bit = [](T field, size_t bit_pos) -> T {
        return (field >> bit_pos) & 1u;
    };

    assert(set_bit(0b0001001, 1, 1) == 0b0001011);
    assert(set_bit(0b0001001, 0, 1) == 0b0001001);
    assert(set_bit(0b0001001, 0, 0) == 0b0001000);
    assert(set_bit(0b0001001, 3, 0) == 0b0000001);
    assert(get_bit(0b0001001, 0) == 1);
    assert(get_bit(0b0001001, 1) == 0);
    assert(get_bit(0b0001001, 2) == 0);
    assert(get_bit(0b0001001, 3) == 1);

    constexpr let set_nth_bit = [](Slice<T>* num, size_t i, T val = 1){
        const size_t digit_i = i / BIT_SIZE<T>;
        const size_t bit_i = i % BIT_SIZE<T>;

        assert(digit_i < num->size);

        T* digit = &(*num)[digit_i];
        *digit = set_bit(*digit, bit_i, val);
    };

    constexpr let get_nth_bit = [](Slice<const T> num, size_t i) -> T {
        const size_t digit_i = i / BIT_SIZE<T>;
        const size_t bit_i = i % BIT_SIZE<T>;

        return get_bit(num[digit_i], bit_i);
    };

    //I am tired of thinking... I am just gonna copy some algorithm from wikipedia
    Slice<T> curr_remainder = trim(*remainder, 0);
    for(size_t i = BIT_SIZE<T> * num.size; i-- > 0; ) 
    {
        let shift_res = shift_up_overflow_batch<T>(&curr_remainder, curr_remainder, 1);
        const T num_ith_bit = get_nth_bit(num, i);

        //this is very disgusting but I am tired
        if(shift_res.overflow != 0 || (num_ith_bit == 1 && curr_remainder.size == 0))
        {
            curr_remainder =  unsafe_slice(curr_remainder, 0, curr_remainder.size + 1);
            *last(&curr_remainder) = shift_res.overflow;

            assert(curr_remainder.size <= remainder->size && "should not go over original remainder boundaries");
        }
        set_nth_bit(remainder, 0, num_ith_bit);

        if(compare<T>(curr_remainder, den) >= 0)
        {
            let sub_res = sub_overflow_batch<T, true>(&curr_remainder, curr_remainder, den);
            assert(sub_res.overflow == 0 && "should not overflow");
            curr_remainder = striped_trailing_zeros<T>(curr_remainder); //for faster checks
                
            if constexpr(DO_QUOTIENT)
                set_nth_bit(quotient, i);
        }
    }

    let stripped_quotient = DO_QUOTIENT 
        ? striped_trailing_zeros<T>(*quotient)
        : trim(*quotient, 0);

    return wrap(Div_Res<T>{stripped_quotient, curr_remainder});
}

template <typename T, 
    bool DO_SHIFT = DO_DEFAULT_DIV_SHIFT_OPTIMS, 
    bool DO_CONST = DO_DEFAULT_DIV_CONST_OPTIMS>
proc rem_bit_by_bit(Slice<T>* remainder, Slice<const T> num, Slice<const T> den) -> Trivial_Maybe<Slice<T>>
{
    static_assert(std::is_unsigned_v<T>);
    Slice<T> quotient = {nullptr, 0};
    let div_res = div_bit_by_bit<T, DO_SHIFT, DO_CONST, false>(&quotient, remainder, num, den);
    if(div_res.has == false)
        return {};

    assert(div_res.quotient.size == 0 && "in case of only remainder the quotient size should be 0");
    assert(div_res.quotient.data == nullptr && "and data nullptr as we set it");
    return wrap(div_res.remainder);
}

//performs to = added + multiplied * coef
//During fuzed quadratic multiply will be called to add back to itself so:
// to += multiplied * coef
template <typename T, bool IN_PLACE = false, bool DO_CONSTS = DO_DEFAULT_MUL_CONST_OPTIMS, bool DO_SHIFT = DO_DEFAULT_MUL_SHIFT_OPTIMS>
proc fuzed_mul_add_overflow_batch(Slice<T>* to, Slice<const T> added, Slice<const T> multiplied, T coeficient,T add_carry = 0, T mul_carry = 0) -> Batch_Overflow<T>
{   
    assert(added.size >= multiplied.size);
    assert(to->size >= added.size);
    assert(check_are_one_way_aliasing<T>(added, *to) == false);
    assert(check_are_one_way_aliasing<T>(multiplied, *to) == false);

    if constexpr(DO_CONSTS)
    {
        const Slice<T> trimmed_to = trim(*to, added.size);
        if(coeficient == 0) 
        {
            if constexpr (IN_PLACE)
                return Batch_Overflow<T>{trim(*to, 0), 0};
            else
            {
                copy_n(to->data, added.data, added.size, Iter_Direction::FORWARD);
                return Batch_Overflow<T>{trimmed_to, 0};
            }
        }

        if(coeficient == 1)
            return add_overflow_batch<T, IN_PLACE>(to, added, multiplied);
    }

    Slice<T> trimmed_to = trim(*to, added.size);
    if(DO_SHIFT && is_power_of_two(coeficient))
    {
        let shift_bits = find_last_set_bit(coeficient);
        const Slice<T> trimmed_to = trim(*to, added.size);
        for (size_t j = 0; j < multiplied.size; j++)
        {
            T curr_multiplied = multiplied[j];
            T curr_add_tp = added[j];

            let mul_res = shift_up_overflow<T>(curr_multiplied, shift_bits, mul_carry);
            let add_res = add_overflow<T>(curr_add_tp, mul_res.value, add_carry);

            trimmed_to[j] = add_res.value;

            mul_carry = mul_res.overflow;
            add_carry = add_res.overflow;
        }
    }
    else
    {
        for (size_t j = 0; j < multiplied.size; j++)
        {
            T curr_multiplied = multiplied[j];
            T curr_add_tp = added[j];

            let mul_res = mul_overflow<T>(curr_multiplied, coeficient, mul_carry);
            let add_res = add_overflow<T>(curr_add_tp, mul_res.value, add_carry);

            trimmed_to[j] = add_res.value;

            mul_carry = mul_res.overflow;
            add_carry = add_res.overflow;
        }
    }
    
    const T combined_carry = add_no_overflow(mul_carry, add_carry);
    return add_overflow_batch_patch<T, IN_PLACE>(&trimmed_to, added, combined_carry, multiplied.size);
}

template <typename T,
    bool DO_SHIFT = DO_DEFAULT_MUL_SHIFT_OPTIMS, 
    bool DO_CONSTS = DO_DEFAULT_MUL_CONST_OPTIMS, 
    bool DO_HALVES = DO_DEFAULT_MUL_HALF_BITS_OPTIMS>
proc mul_quadratic(Slice<T>* to, Slice<T>* temp, Slice<const T> left, Slice<const T> right) -> void
{
    static_assert(std::is_unsigned_v<T>);
    assert(check_are_aliasing<T>(*to, left) == false);
    assert(check_are_aliasing<T>(*to, right) == false);
    assert(check_are_aliasing<T>(*temp, left) == false);
    assert(check_are_aliasing<T>(*temp, right) == false);
    assert(check_are_aliasing<T>(*to, *temp) == false);

    if(left.size < right.size)
        std::swap(left, right);

    const size_t max_result_size = left.size + right.size + 1;
    assert(to->size >= max_result_size);
    assert(temp->size >= right.size + 1);
    null_n(to->data, to->size);

    for(size_t i = 0; i < right.size; i++)
    {
        let mul_res = mul_overflow_batch<T, DO_SHIFT, DO_CONSTS, DO_HALVES>(temp, left, right[i]);

        let temp_num = unwrap(to_number<u64>(mul_res.slice));

        //mul_overflow_batch is limited to left.size elems => in case of overflow
        // even if the overflow would fit into temp it is not added => we add it back
        let mul_size = mul_res.slice.size;
        let mul_slice = trim(*temp, mul_size + 1); 
        mul_slice[mul_size] = mul_res.overflow;

        mut shifted_to = slice(*to, i);
        assert(shifted_to.size >= mul_slice.size);

        let add_res = add_overflow_batch<T, true>(&shifted_to, shifted_to, mul_slice);
        assert(add_res.overflow == 0 && "no final carry should be left");
    }
}

template <typename T, bool DO_CONSTS = DO_DEFAULT_MUL_CONST_OPTIMS>
proc mul_quadratic_fused(Slice<T>* to, Slice<const T> left, Slice<const T> right) -> void
{
    static_assert(std::is_unsigned_v<T>);
    assert(check_are_aliasing<T>(*to, left) == false);
    assert(check_are_aliasing<T>(*to, right) == false);

    if(left.size < right.size)
        std::swap(left, right);

    const size_t max_result_size = left.size + right.size + 1;
    assert(to->size >= max_result_size);
    null_n(to->data, to->size);

    for(size_t i = 0; i < right.size; i++)
    {
        mut shifted_to = slice(*to, i);
        const T curr_right = right[i];

        let mul_add_res = fuzed_mul_add_overflow_batch<T, true, DO_CONSTS>(&shifted_to, shifted_to, left, curr_right);
        assert(mul_add_res.overflow == 0 && "no carry should happen in this case");
    }
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
template <typename From>
runtime_func required_size_to_base(size_t num_size, size_t to_base) -> size_t
{
    double in_one_digit = round(pow(2, BIT_SIZE<From>));
    return cast(size_t) ceil(size_to_fit_base(in_one_digit, cast(double) num_size, cast(double) to_base));
}

//what is the size of the output buffer when converting from base rep to number?
template <typename From>
runtime_func required_size_from_base(size_t represenation_size, size_t from_base) -> size_t
{
    double in_one_digit = round(pow(2, BIT_SIZE<From>));
    return cast(size_t) ceil(size_to_fit_base(cast(double) from_base, cast(double) represenation_size, in_one_digit));
}

template <typename From, typename To>
proc to_base(Slice<To>* rep, Slice<From>* temp, Slice<const From> num, To base) -> Slice<To>
{
    //init
    static_assert(std::is_unsigned_v<From>);
    //static_assert(sizeof(To) <= sizeof(From));
    if constexpr(std::is_same_v<From, To>)
        assert(check_are_aliasing<From>(num, *rep) == false);

    const From cast_base = cast(From) base;

    assert(is_striped_form(num));
    assert(temp->size >= num.size && "temp must be big enough");
    assert(base >= 2 && "base must be bigger than two");
    assert(cast(To) cast_base == base && "base must fit in type");
    assert(high_bits(cast_base) == 0 && "base must be low bits only of From type (so that short div algorhimt can be used)");

    if constexpr (DO_RUNTIME_ONLY_CHECKS)
    {
        size_t at_least = required_size_to_base<From>(num.size, base);
        assert(at_least <= rep->size && "there must be enough size rep represent the number even in the worst case scenario");
    }

    size_t to_size = 0;
    if(num.size == 0)
        return trim(*rep, to_size);

    //convert
    Slice<From> val_left = *temp;
    while(true)
    {
        const Slice<const From> curr_div_from = to_size == 0 
            ? num : val_left;
            
        if(curr_div_from.size == 0)
            break;

        let res = div_overflow_low_batch<From>(&val_left, curr_div_from, cast_base);
        val_left = striped_trailing_zeros(val_left);
        let digit = cast(To) res.overflow;
        (*rep)[to_size] = digit;
        to_size++;
    }


    //finish
    Slice<To> converted = trim(*rep, to_size);
    assert(is_striped_form(converted));
    reverse(&converted);
    return converted;
}

template <typename From, typename To>
proc from_base(Slice<To>* num, Slice<const From> rep, To base) -> Slice<To>
{
    static_assert(std::is_unsigned_v<To>);
    assert(is_reverse_striped_form(rep));
    if constexpr(std::is_same_v<From, To>)
        assert(check_are_aliasing<From>(rep, *num) == false);

    if constexpr (DO_RUNTIME_ONLY_CHECKS)
    {
        size_t at_least = required_size_from_base<From>(rep.size, base);
        assert(at_least <= num->size && "there must be enough size num represent the number even in the worst case scenario");
    }

    constexpr size_t single_num_digits = digits_to_represent<To, From>();
    To single_num_vals[single_num_digits] = {0};
    Slice<To> single_num = {single_num_vals, single_num_digits};

    null_n(num->data, num->size);

    for(size_t i = 0; i < rep.size; i++)
    {
        let digit = rep[i];
        let converted = from_number<To, From>(&single_num, digit);
        let mul_res = mul_overflow_batch<To>(num, *num, base);
        let add_res = add_overflow_batch<To, true>(num, *num, converted);

        assert(add_res.overflow == 0);
        assert(mul_res.overflow == 0);
    }

    return striped_trailing_zeros(*num);
}