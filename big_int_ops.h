#pragma once

#include "preface.h"

using Max_Unsigned_Type = u64;

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
constexpr size_t BIT_SIZE = sizeof(T) * CHAR_BIT;
template <typename T>
constexpr size_t HALF_BIT_SIZE = BIT_SIZE<T> / 2;

template <integral T>
constexpr T FULL_MASK = cast(T) Max_Unsigned_Type(-1);
/*
template <typename T>
struct Span
{
    T* data;
    size_t size;

    func& operator[](size_t index) const noexcept { assert(index < size && "index out of range"); return data[index]; }
    func& operator[](size_t index) noexcept       { assert(index < size && "index out of range"); return data[index]; }
};
*/

//@TODO: Tidy up the div proc - move bit manipulation out
//@TODO: Merge typed and utyped tests so that the typed part executes only when the supplied typed matches
//@TODO: Figure out the smallest possible requirements for div
//@TODO: Add alias checking function and alias check functions that donst support self assign
//@TODO: Add asserts for aliasing that is either full or none
//@TODO: Instantiate throwing versions of functions as well
//@TODO: Insatntiate different optims
//@TODO: Add more fine grained optim switches
//@TODO: Create interface for parsing and printing functions (taking functions probably or returnign state to be executed repeatedly) 
//@TODO: Hook to the existing parsing functions I have made
//@TODO: Make Better multiply algorhitm
//@TODO: Maybe just copy over 
//@TODO: Add define flags setting default
//@TODO: Remove extra dependencies
//@TODO: Make c++ 17 friendly (remove concepts add type asserts)
//@TODO: Figure out param passing

static constexpr bool DO_MUL_OPTIMS = true;
static constexpr bool DO_DIV_OPTIMS = true;
static constexpr bool DO_DEFAULT_THROW = true;

template <integral T>
func high_mask(size_t index = HALF_BIT_SIZE<T>) -> T {
    return cast(T) (FULL_MASK<T> << index);
}

template <integral T>
func low_mask(size_t index = HALF_BIT_SIZE<T>) -> T {
    return cast(T) ~high_mask<T>(index);
}

template <integral T>
func high_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
    return value >> index;
}

template <integral T>
func low_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
    return value & low_mask<T>(index);
};  

template <integral T>
func combine_bits(T low, T high, size_t index = HALF_BIT_SIZE<T>) -> T {
    return low_bits(low, index) | (high << index);
};

template <integral T>
func dirty_combine_bits(T low, T high, size_t index = HALF_BIT_SIZE<T>) -> T {
    assert(high_bits(low, index) == 0 && "low must not have high bits use combine_bits instead");
    return low | (high << index);
};

template <integral T>
func integral_log2(T val) -> T  
{
    if constexpr(sizeof(T) > 8)
    {
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

template <integral T>
struct Overflow
{
    T value;
    T overflow;

    bool constexpr operator ==(Overflow const&) const noexcept = default;
};

template <integral T>
struct Batch_Overflow
{
    size_t size; //count of digits written to output span
    T overflow;

    bool constexpr operator ==(Batch_Overflow const&) const noexcept = default;
};


template <integral T>
func zeros_from_index(span<T> num) -> size_t
{
    size_t to = num.size();
    for(; to-- > 0;)
        if(num[to] != 0)
            break;

    return to + 1;
}

template <integral T>
func striped_trailing_zeros(span<T> num) -> span<T>
{
    return num.subspan(0, zeros_from_index(num));
}

template <integral T>
func is_striped_form(span<T> num) -> bool
{
    return zeros_from_index(num) == num.size();
}

template <integral Digit, integral Number>
func digits_to_represent() -> size_t
{
    return (sizeof(Digit) + sizeof(Number) - 1) / sizeof(Number);
}

template <integral To = Max_Unsigned_Type, integral T = u8>
func to_number(span<T> bignum) -> Trivial_Maybe<To>
{
    let stripped = striped_trailing_zeros(bignum);

    const size_t total_size = stripped.size() * sizeof(T);
    if(total_size > sizeof(To))
        return Trivial_Maybe<To>{false};

    To out = 0;
    for(size_t i = 0; i < stripped.size(); i++)
        out |= cast(To)(stripped[i]) << (i * BIT_SIZE<T>);

    return wrap<To>(out);
}

template <integral T, integral From = Max_Unsigned_Type>
proc from_number(span<T>* bignum, From from) -> size_t
{
    constexpr size_t count = digits_to_represent<T, From>();
    assert(bignum.size() >= count);

    size_t i = 0;
    for(; i < count; i++)
    {
        T curr_digit = from >> (i * BIT_SIZE<T>);
        if(curr_digit == 0)
            break;

        (*bignum)[i] = cast(T) curr_digit;
    }

    return i;
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


template <integral T>
func add_overflow(T left, T right, T carry = 0) -> Overflow<T>
{
    assert(high_bits(carry) == 0 && "carry must be half bits use add_overflow_any instead");
    const T res_low = low_bits(left) + low_bits(right) + carry;
    const T middle_carry = high_bits(res_low);
    
    const T res_high = high_bits(left) + high_bits(right) + middle_carry;
    const T new_carry = high_bits(res_high);

    return Overflow<T>{combine_bits(res_low, res_high), new_carry};
}

namespace detail
{
    template <integral T, integral... Ts>
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
        constexpr size_t max_args = cast(size_t) 1 << HALF_BIT_SIZE<T>;
        static_assert(arg_count <= max_args && "too many operands");

        const T res_low = (low_bits(cast(T) args) + ...);
        const T middle_carry = high_bits(res_low);

        const T res_high = (high_bits(cast(T) args) + ...) + middle_carry;
        const T new_carry = high_bits(res_high);

        return Overflow<T>{combine_bits(res_low, res_high), new_carry};
    }

    template <integral T, integral... Ts>
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

template <integral T, integral... Ts>
func add_overflow_any(T left, T right, Ts... rest) -> Overflow<T>
{
    //we separate it like so that the user must supply at least 2 operands
    return detail::add_overflow_any<T>(left, right, rest...);
}

template <integral T, integral... Ts>
func add_no_overflow(T left, T right, Ts... rest) -> T
{
    return detail::add_no_overflow<T>(left, right, rest...);
}

template <integral T>
func sub_overflow(T left, T right, T carry = 0) -> Overflow<T>
{
    const T res_low = low_bits(left) - low_bits(right) - carry;
    const T middle_carry = res_low >> (BIT_SIZE<T> - 1);

    const T res_high = high_bits(left) - high_bits(right) - middle_carry;
    const T new_carry = res_high >> (BIT_SIZE<T> - 1);

    return Overflow<T>{combine_bits(res_low, res_high), new_carry};
}

template <integral T>
func add_carry_patch_step(T left, T carry) -> Overflow<T>
{
    T res = left + carry;
    return Overflow<T>{res, cast(T)(res < left)};
}

template <integral T>
func sub_carry_patch_step(T left, T carry) -> Overflow<T>
{
    T res = left - carry;
    return Overflow<T>{res, cast(T)(res > left)};
}

template <integral T>
func add_overflow_batch(span<T>* to, span<const T> left, span<const T> right, T carry_in = 0) -> Batch_Overflow<T>
{
    assert(to->size() >= right.size());
    assert(to->size() >= left.size());
    assert(high_bits(carry_in) == 0);

    if(left.size() < right.size())
        std::swap(left, right);

    T carry = carry_in;
    proc update = [&](size_t i, Overflow<T> oveflow) {
        (*to)[i] = oveflow.value;
        carry = oveflow.overflow;
    };

    const size_t min_size = min(left.size(), right.size());
    const size_t max_size = max(left.size(), right.size());
    for (size_t i = 0; i < min_size; i++)
        update(i, add_overflow<T>(left[i], right[i], carry));

    //@TODO: add early stop optim
    for (size_t i = right.size(); i < left.size(); i++)
        update(i, add_carry_patch_step<T>(left[i], carry));

    return Batch_Overflow<T>{max_size, carry};
}

template <integral T>
func sub_overflow_batch(span<T>* to, span<const T> left, span<const T> right, T carry_in = 0) -> Batch_Overflow<T>
{
    assert(to->size() >= right.size());
    assert(to->size() >= left.size());
    assert(high_bits(carry_in) == 0);

    T carry = carry_in;
    proc update = [&](size_t i, Overflow<T> oveflow) {
        (*to)[i] = oveflow.value;
        carry = oveflow.overflow;
    };

    const size_t min_size = min(left.size(), right.size());
    const size_t max_size = max(left.size(), right.size());
    for (size_t i = 0; i < min_size; i++)
        update(i, sub_overflow<T>(left[i], right[i], carry));

    for (size_t i = left.size(); i < right.size(); i++)
        update(i, sub_overflow<T>(0, right[i], carry));

    for (size_t i = right.size(); i < left.size(); i++)
        update(i, sub_carry_patch_step<T>(left[i], carry));

    return Batch_Overflow<T>{max_size, carry};
}

template <integral T>
func complement_overflow(T val, T carry) -> Overflow<T>
{
    const T val_inv = ~val;

    const T complement_low = low_bits(val_inv) + carry;
    const T middle_carry = high_bits(complement_low);

    const T complement_high = high_bits(val_inv) + middle_carry;
    const T new_carry = high_bits(complement_high);

    return Overflow<T>{combine_bits(complement_low, complement_high), new_carry};
}

template <integral T>
proc complement_overflow_batch(span<T>* to, span<const T> left, T carry_in = 1) -> Batch_Overflow<T>
{
    assert(to->size() >= left.size());
    assert(is_striped_form(left));
    assert(carry_in == 1 || carry_in == 0);

    T carry = carry_in;
    for(size_t i = 0; i < left.size(); i++)
    {
        let res = complement_overflow<T>(left[i], carry);
        (*to)[i] = res.value;
        carry = res.overflow;
    }

    return Batch_Overflow<T>{left.size(), carry};
}
template <integral T>
func shift_up_overflow(T low_item, size_t by_bits, T high_item) -> Overflow<T>
{
    assert(by_bits < BIT_SIZE<T>);
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

    const T low = high_bits(low_item, remaining_bits);
    const T high = low_bits(high_item, remaining_bits);

    const T composed = low | high << by_bits;

    const T out = low_bits(high_item, remaining_bits);
    const T shifted_out = out;
    return Overflow<T>{composed, shifted_out};
}

template <integral T>
func shift_down_overflow(T low_item, size_t by_bits, T high_item) -> Overflow<T>
{
    assert(by_bits < BIT_SIZE<T>);
    const size_t remaining_bits = BIT_SIZE<T> - by_bits;

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

template <integral T>
proc shift_up_overflow_batch(span<T>* out, span<const T> in, size_t by_bits, T carry_in = 0) -> Batch_Overflow<T>
{
    assert(out->size() >= in.size());
    assert(by_bits < BIT_SIZE<T>);

    if(in.size() == 0)
        return Batch_Overflow<T>{in.size(), shift_up_overflow<T>(carry_in, by_bits, 0).value};

    const T shifted_out = in[in.size() - 1];
    T prev = shifted_out;
    for (size_t i = in.size(); i-- > 1; )
    {
        let curr = in[i - 1];
        let res = shift_up_overflow<T>(curr, by_bits, prev);
        (*out)[i] = res.value;
        prev = curr;
    }

    {
        let res = shift_up_overflow<T>(carry_in, by_bits, prev);
        (*out)[0] = res.value;
    }

    let carry_out = shift_up_overflow<T>(shifted_out, by_bits, 0);
    return Batch_Overflow<T>{in.size(), carry_out.value};
}


template <integral T>
proc shift_down_overflow_batch(span<T>* out, span<const T> in, size_t by_bits, T carry_in = 0) -> Batch_Overflow<T>
{
    assert(out->size() >= in.size());
    assert(by_bits < BIT_SIZE<T>);

    if(in.size() == 0)
        return Batch_Overflow<T>{in.size(), shift_down_overflow<T>(0, by_bits, carry_in).value};

    const T shifted_out = in[0];
    T prev = shifted_out;
    for (size_t i = 0; i < in.size() - 1; i++)
    {
        let curr = in[i + 1];
        let res = shift_down_overflow<T>(prev, by_bits, curr);
        (*out)[i] = res.value;
        prev = curr;
    }

    {
        let res = shift_down_overflow<T>(prev, by_bits, carry_in);
        (*out)[in.size() - 1] = res.value;
    }

    let carry_out = shift_down_overflow<T>(0, by_bits, shifted_out);
    return Batch_Overflow<T>{in.size(), carry_out.value};
}

//return -1 if left < right
//        0 if left == right
//        1 if left > right   
template <integral T>
func compare(span<const T> left, span<const T> right) -> int
{
    assert(is_striped_form(left));
    assert(is_striped_form(right));
    if(left.size() < right.size())
        return -1;

    if(left.size() > right.size())
        return 1;

    for(size_t i = left.size(); i-- > 0;)
    {
        if(left[i] < right[i])
            return -1;
        if(left[i] > right[i])
            return 1;
    }

    return 0;
}

enum class Mul_Overflow_Optims
{
    NONE,
    HIGH_BITS_ONLY,
    LOW_BITS_ONLY,
};

template <integral T, Mul_Overflow_Optims OPTIMS = Mul_Overflow_Optims::NONE>
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
    T low_mixed = low_bits(s_mixed_ns) << h; // even though we are taking lower half of s_mixed we have to insert it to upper half
    T high_mixed = high_bits(s_mixed_ns) + (s_mixed_ns_overflow << h); //add overflow as another digit

    //distribute mixed to appropriate halves
    let curr_value = add_overflow_any<T>(s_cur, low_mixed, last_value); //also add last_value to save ops
    T next_value =   add_no_overflow<T>(s_next_ns, high_mixed, curr_value.overflow); //also add the overflow

    return Overflow<T>{curr_value.value, next_value};
}

template <integral T>
func is_power_of_two(T num)
{
    return num != 0 && (num & (num - 1)) == 0;
}

template <integral T, bool DO_OPTIMS = DO_MUL_OPTIMS>
proc mul_overflow_batch(span<T>* to, span<const T> left, T right, T carry_in = 0) -> Batch_Overflow<T>
{
    assert(to->size() >= left.size());

    //@TODO: Add more fine grained optim switches
    //@TODO: rework default optim naming
    if(right == 0) 
        return Batch_Overflow<T>{0, 0};

    T carry = carry_in;
    if constexpr(DO_OPTIMS)
    {
        if(right == 1)
        {
            safe_copy_n(to->data(), left.data(), left.size());
            return Batch_Overflow<T>{left.size(), 0};
        }

        if(is_power_of_two(right))
        {
            let shift_bits = integral_log2(right);
            return shift_up_overflow_batch(to, left, shift_bits);
        }

        if(high_bits(right) == 0)
        {
            for (size_t i = 0; i < left.size(); i++)
            {
                let res = mul_overflow<T, Mul_Overflow_Optims::LOW_BITS_ONLY>(left[i], right, carry);
                (*to)[i] = res.value;
                carry = res.overflow;
            }
            return Batch_Overflow<T>{left.size(), carry};
        }

        if(low_bits(right) == 0)
        {
            for (size_t i = 0; i < left.size(); i++)
            {
                let res = mul_overflow<T, Mul_Overflow_Optims::HIGH_BITS_ONLY>(left[i], right, carry);
                (*to)[i] = res.value;
                carry = res.overflow;
            }
            return Batch_Overflow<T>{left.size(), carry};
        }
    }

    for (size_t i = 0; i < left.size(); i++)
    {
        let res = mul_overflow<T>(left[i], right, carry);
        (*to)[i] = res.value;
        carry = res.overflow;
    }

    return Batch_Overflow<T>{left.size(), carry};
}


template <integral T>
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

struct Div_Res
{
    size_t quotient_size;
    size_t remainder_size;
};

namespace detail
{
    template <integral T, bool DO_OPTIMS = DO_DIV_OPTIMS>
    proc div_overflow_low_batch(span<T>* to, span<const T> left, T right, T carry_in = 0) -> Batch_Overflow<T>
    {
        assert(to->size() >= left.size());
        assert(high_bits(right) == 0 && "only works for divisors under half bit size");
        assert(right != 0 && "cannot divide by zero");

        if constexpr(DO_OPTIMS)
        {
            if(right == 1)
            {
                safe_copy_n(to->data(), left.data(), left.size());
                return Batch_Overflow<T>{left.size(), 0}; 
            }

            if(is_power_of_two(right))
            {
                let shift_bits = integral_log2(right);
                let shift_res = shift_down_overflow_batch(to, left, shift_bits);
                //because of the way shifting down result is defined we have to process the reuslt
                const T carry = shift_res.overflow >> (BIT_SIZE<T> - shift_bits);
                return Batch_Overflow<T>{left.size(), carry};
            }
        }

        T carry = carry_in;
        for (size_t i = left.size(); i-- > 0; )
        {
            let res = div_overflow_low<T>(left[i], right, carry);
            (*to)[i] = res.value;
            carry = res.overflow;
        }

        return Batch_Overflow<T>{left.size(), carry};

    }

    template <integral T, bool DO_OPTIMS = DO_DIV_OPTIMS, bool DO_QUOTIENT = true>
    proc div_bit_by_bit(span<T>* quotient, span<T>* remainder, span<const T> num, span<const T> den) -> Div_Res
    {
        assert(is_striped_form(num));
        assert(is_striped_form(den));
        assert(den.size() != 0 && "cannot divide by zero");

        const size_t min_size = min(num.size(), den.size());
        assert(remainder->size() >= min_size);
        assert(quotient->size() >= min_size);

        null_n(quotient->data(), quotient->size());
        null_n(remainder->data(), remainder->size());

        if(num.size() < den.size())
        {
            safe_copy_n(remainder->data(), num.data(), num.size());
            return {0, num.size()};
        }

        if constexpr(DO_OPTIMS)
        {
            if(den.size() == 1 && high_bits(den.back()) == 0)
            {
                let div_res = div_overflow_low_batch<T, DO_OPTIMS>(quotient, num, den.back());
                assert(remainder->size() > 0 && "at this point should be 0");

                size_t remainder_size = 0;
                if(div_res.overflow != 0)
                {
                    (*remainder)[0] = div_res.overflow;
                    remainder_size = 1;
                }

                let stripped_quotient = striped_trailing_zeros<T>(*quotient);
                return {stripped_quotient.size(), remainder_size};
            }
        }

        constexpr let set_bit = [](T field, size_t bit_pos, T val) -> T {
            assert(val == 0 || val == 1);

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

        constexpr let set_nth_bit = [](span<T>* num, size_t i, T val = 1){
            const size_t digit_i = i / BIT_SIZE<T>;
            const size_t bit_i = i % BIT_SIZE<T>;

            assert(digit_i < num->size());

            T& digit = (*num)[digit_i];
            digit = set_bit(digit, bit_i, val);
        };

        constexpr let get_nth_bit = [](span<const T> num, size_t i) -> T {
            const size_t digit_i = i / BIT_SIZE<T>;
            const size_t bit_i = i % BIT_SIZE<T>;

            return get_bit(num[digit_i], bit_i);
        };

        //I am tired of thinking... I am just gonna copy some algorithm from wikipedia
        span<T> curr_remainder = span<T>(remainder->data(), 0);
        for(size_t i = BIT_SIZE<T> * num.size(); i-- > 0; ) 
        {
            let shift_res = shift_up_overflow_batch<T>(&curr_remainder, curr_remainder, 1);
            const T num_ith_bit = get_nth_bit(num, i);

            //this is very disgusting but I am tired
            if(shift_res.overflow != 0 || (num_ith_bit == 1 && curr_remainder.size() == 0))
            {
                curr_remainder = span<T>(curr_remainder.data(), curr_remainder.size() + 1);
                curr_remainder.back() = shift_res.overflow;

                assert(curr_remainder.size() <= remainder->size() && "should not go over original remainder boundaries");
            }
            set_nth_bit(remainder, 0, num_ith_bit);

            if(compare<T>(curr_remainder, den) >= 0)
            {
                let sub_res = sub_overflow_batch<T>(&curr_remainder, curr_remainder, den);
                assert(sub_res.overflow == 0 && "should not overflow");
                curr_remainder = striped_trailing_zeros<T>(curr_remainder); //for faster checks
                
                //in case we only want modulo
                if constexpr(DO_QUOTIENT)
                    set_nth_bit(quotient, i);
            }
        }


        let stripped_quotient = striped_trailing_zeros<T>(*quotient);
        return {stripped_quotient.size(), curr_remainder.size()};
    }

}

struct Divide_By_Zero_Exception{};

template <integral T, bool DO_THROW = DO_DEFAULT_THROW, bool DO_OPTIMS = DO_DIV_OPTIMS>
proc div_overflow_low_batch(span<T>* to, span<const T> left, T right, T carry_in = 0) -> std::conditional_t<DO_THROW, Batch_Overflow<T>, Trivial_Maybe<Batch_Overflow<T>>>
{
    using Maybe = Trivial_Maybe<Batch_Overflow<T>>;

    if(right == 0)
    {
        if constexpr(DO_THROW)
            throw Divide_By_Zero_Exception{};
        else
            return Maybe{false};
    }

    let res = detail::div_overflow_low_batch<T, DO_OPTIMS>(to, left, right, carry_in);
    if constexpr(DO_THROW)
        return res;
    else
        return Maybe{true, res};
}

template <integral T, bool DO_THROW = DO_DEFAULT_THROW, bool DO_OPTIMS = DO_DIV_OPTIMS>
proc div_bit_by_bit(span<T>* quotient, span<T>* remainder, span<const T> left, span<const T> right) -> std::conditional_t<DO_THROW, Div_Res, Trivial_Maybe<Div_Res>>
{
    using Maybe = Trivial_Maybe<Div_Res>;
    assert(is_striped_form(right));
    if(right.size() == 0)
    {
        if constexpr(DO_THROW)
            throw Divide_By_Zero_Exception{};
        else
            return Maybe{false};
    }

    let res = detail::div_bit_by_bit<T, DO_OPTIMS>(quotient, remainder, left, right);
    if constexpr(DO_THROW)
        return res;
    else
        return Maybe{true, res};
}

template <integral T>
proc mul_quadratic(span<T>* to, span<T>* temp, span<const T> left, span<const T> right) -> void
{
    if(left.size() < right.size())
        std::swap(left, right);

    size_t max_result_size = left.size() + right.size() + 1;
    assert(to->size() >= max_result_size);
    assert(temp->size() >= right.size() + 1);
    null_n(to->data(), to->size());

    for(size_t i = 0; i < right.size(); i++)
    {
        let mul_res = mul_overflow_batch<T>(temp, left, right[i]);
        size_t curr_size = mul_res.size + 1;
        (*temp)[mul_res.size] = mul_res.overflow;
        
        mut shifted_to = to->subspan(i);
        assert(shifted_to.size() >= curr_size);
        assert(temp->size() >= curr_size);

        let add_res = add_overflow_batch<T>(&shifted_to, shifted_to, temp->subspan(0, curr_size));
        assert(add_res.overflow == 0); //no carry should happen in this case
    }
}


template <integral T>
proc mul_quadratic_fused(span<T>* to, span<const T> left, span<const T> right) -> void
{
    if(left.size() < right.size())
        std::swap(left, right);

    size_t max_result_size = left.size() + right.size() + 1;
    assert(to->size() >= max_result_size);
    null_n(to->data(), to->size());

    for(size_t i = 0; i < right.size(); i++)
    {
        mut shifted_to = to->subspan(i);
        const T curr_right = right[i];
        if(curr_right == 0) 
            continue;

        if constexpr(DO_MUL_OPTIMS)
        {
            if(curr_right == 1)
            {
                let add_res = add_overflow_batch<T>(&shifted_to, span<const T>{shifted_to}, left);
                assert(add_res.overflow == 0); //no carry should happen in this case
                continue;
            }
        }

        T mul_carry = 0;
        T add_carry = 0;
        for (size_t j = 0; j < left.size(); j++)
        {
            T curr_left = left[j];
            T curr_to = shifted_to[j];

            let mul_res = mul_overflow<T>(curr_left, curr_right, mul_carry);
            let add_res = add_overflow<T>(curr_to, mul_res.value, add_carry);

            shifted_to[j] = add_res.value;

            mul_carry = mul_res.overflow;
            add_carry = add_res.overflow;
        }

        const T combined_carry = add_no_overflow(mul_carry, add_carry);
        if(combined_carry != 0)
        {
            for (size_t j = left.size(); j < shifted_to.size(); j++)
            {
                let patch_res = add_carry_patch_step<T>(shifted_to[j], combined_carry);
                shifted_to[j] = patch_res.value;
                if(patch_res.overflow == 0)
                    break;
            }
        }
    }
}

//template <integral From, integral To>
//proc to_base(span<const From> num, span<From*> div_temp1, span<From*> div_temp2, span<From*> remainder_temp, span<To>* to, size_t base, let conversion_func)
//{
//    constexpr size_t item_count = digits_to_represent<From, size_t>();
//    
//    From base_rep_arr[item_count];
//    span<From> whole_base_rep = base_rep_arr;
//
//    const size_t required_size = from_number(&whole_base_rep, base);
//    const span<From> base_rep = span<From>(base_rep_arr, required_size)
//
//
//
//    //10541 -> 1054 | 1
//}
