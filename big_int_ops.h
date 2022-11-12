#pragma once

#include <algorithm>
#include <ranges>
#include "big_int.h"

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

template <integral T>
func high_mask(size_t index = HALF_BIT_SIZE<T>) -> T
{
    return cast(T) (FULL_MASK<T> << index);
}

template <integral T>
func low_mask(size_t index = HALF_BIT_SIZE<T>) -> T
{
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
struct Overflow
{
    T value;
    T overflow;
};

template <integral T>
struct Batch_Overflow
{
    size_t size; //count of digits written to output span
    T overflow;
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

template <typename T>
func as_number(span<T> nums) -> Max_Unsigned_Type
{
    constexpr size_t rep_size = sizeof(Max_Unsigned_Type) / sizeof(T);

    let stripped = striped_trailing_zeros(nums);

    const size_t total_size = stripped.size() * sizeof(T);
    assert(total_size <= sizeof(Max_Unsigned_Type));

    Max_Unsigned_Type out = 0;
    for(size_t i = 0; i < stripped.size(); i++)
        out |= cast(Max_Unsigned_Type)(stripped[i]) << (i * sizeof(T) * CHAR_BIT);

    return out;
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
//  (and we would also get about 2x - 4x speedup)


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
    assert(carry_in == 1 || carry_in == 0);

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
    assert(carry_in == 1 || carry_in == 0);

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

/*
//shifts up by `by_bits` the `in` digits handling carry between the number. Can shift up by any number of bits even surpassing the
// bit size of individual digits. Accepts carry_in which is shifted into the first (smallest) digit not completely shifted out. 
// (ie when shifting by 9 span of u8's the shift can be though of as shifting one digit forward and then doing bitshift by 1 on the whole number again.
//  The carry_in is used only for that bitshift not for digit shift)
// Returns the last digits shifted out part as a carry (again with the same logic while shifting more than one digits bit size as with carry_in)
template <integral T>
func shift_up_overflow_batch(span<T>* out, span<const T> in, size_t by_bits, size_t out_bounds, T carry_in) -> Batch_Overflow<T>
{
    assert(out->size() >= out_bounds);
    assert(out->size() >= in.size());
    size_t shift_items_low = by_bits / BIT_SIZE<T>; 
    size_t shift_items_high = div_round_up(by_bits, BIT_SIZE<T>);
    size_t shift_bits = by_bits - shift_items_low * BIT_SIZE<T>;
    size_t remaining_bits = BIT_SIZE<T> - shift_bits;
    size_t to_index = min(out_bounds, in.size() + shift_items_high);

    T high_item = 0;
    T last_shifted = 0;
    if(to_index >= shift_items_high + 1)
    {
        last_shifted = in[to_index - shift_items_high - 1];

        //when we limit the written to index by out_bounds the first actually written digit still needs
        // to have its upper half talken from somewhere. We in that case set it to the following item which happens
        // to be the last_shifted. Again this only happens when the size is bound by out_bounds
        if(to_index == out_bounds)
            high_item = last_shifted;
    }

    for (size_t i = to_index; i-- > shift_items_high; )
    {
        T low_item = in[i - shift_items_high];

        T high = low_bits(high_item, remaining_bits);
        T low = high_bits(low_item, remaining_bits);

        T composed = low | high << shift_bits;

        (*out)[i] = composed;
        high_item = low_item;
    }

    //we patch the first digit individually by 'shifting in' the carry_in 
    {
        T high = low_bits(high_item, remaining_bits);
        T low = high_bits(carry_in, remaining_bits);

        T composed = low | high << shift_bits;
        (*out)[shift_items_low] = composed;
    }

    null_n(out->data(), min(out_bounds, shift_items_low));

    //modulo for the rare case that we shift exactly by digits bit size
    T carry_out = high_bits(last_shifted, remaining_bits % BIT_SIZE<T>);

    return Batch_Overflow<T>{carry_out, out_bounds};
}
*/

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

    int diff = cast(int) left.size() - cast(int) right.size();
    if(diff != 0)
        return diff;

    for(size_t i = 0; i < left.size(); i++)
    {
        if(left[i] < right[i])
            return -1;
        if(left[i] > right[i])
            return 1;
    }

    return 0;
}


template <integral T>
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
    T s_cur = l1*r1;
    let s_mixed_ns = add_overflow<T>(l1*r2, l2*r1);
    T s_next_ns = l2*r2;

    //split mixed
    T low_mixed = low_bits(s_mixed_ns.value) << h; // even though we are taking lower half of s_mixed we have to insert it to upper half
    T high_mixed = high_bits(s_mixed_ns.value) + (s_mixed_ns.overflow << h); //add overflow as another digit

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

static constexpr bool DO_MUL_OPTIMS = true;
template <integral T>
proc mul_overflow_batch(span<T>* to, span<const T> left, T right, T carry_in = 0) -> Batch_Overflow<T>
{
    assert(to->size() >= left.size());

    if(right == 0) 
        return Batch_Overflow<T>{0, 0};

    if constexpr(DO_MUL_OPTIMS)
    {
        if(right == 1)
        {
            safe_copy_n(to->data(), left.data(), left.size());
            return Batch_Overflow<T>{left.size(), 0};
        }

        /*if(is_power_of_two(right))
        {
            let res = shift_up_overflow_batch(to, left, right);
            return Batch_Overflow<T>{left.size(), res};
        }*/
    }

    T last_value = carry_in;
    for (size_t i = 0; i < left.size(); i++)
    {
        let res = mul_overflow<T>(left[i], right, last_value);
        (*to)[i] = res.value;
        last_value = res.overflow;
    }

    let obtained = as_number(*to);

    return Batch_Overflow<T>{
        left.size(),
        last_value,
    };
}

template <integral T>
proc mul_quadratic(span<T>* to, span<T>* temp, span<const T> left, span<const T> right) -> void
{
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

template <integral From, integral To>
proc to_base(span<const From> num, span<To>* to, size_t base, let conversion_func)
{
}
