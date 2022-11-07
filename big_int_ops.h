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
template <integral T>
constexpr T LOW_MASK = cast(T)(FULL_MASK<T> >> HALF_BIT_SIZE<T>);
template <integral T>
constexpr T HIGH_MASK = cast(T)(FULL_MASK<T> << HALF_BIT_SIZE<T>);
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
template <typename T>
func as_number(span<T> nums) -> Max_Unsigned_Type
{
    constexpr size_t rep_size = sizeof(Max_Unsigned_Type) / sizeof(T);

    size_t total_size = nums.size() * sizeof(T);
    assert(total_size <= sizeof(Max_Unsigned_Type));

    Max_Unsigned_Type out = 0;
    for(size_t i = 0; i < nums.size(); i++)
        out |= cast(Max_Unsigned_Type)(nums[i]) << (i * sizeof(T) * CHAR_BIT);

    return out;
}

//TODO: remove dependance on std::span instead provide a trivially constructibe span struct - the goal of this is to move the efective core of the library away from any particular
//       implementation.
//      Add a query option for every function that will return the required min size for the operation
//      Maybe wrap all functions inside a struct with type T so we dont have to template evrything => easier testing => harder to extend harder to import into namespace
//      Remove integral constrait on functions so we can compile with c++17 (add static asserts instead)

//0xFFAABB, 16 => 0xAABB => >> index
template <integral T>
func high_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
    return value >> index;
}

template <integral T>
func low_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
    return value & (FULL_MASK<T> >> index);
};  

template <integral T>
func combine_bits(T low, T high, size_t index = HALF_BIT_SIZE<T>) -> T {
    return (low & (FULL_MASK<T> >> index)) | (high << index);
};

template <integral T>
struct Carried
{
    T value;
    T carry;
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

template <integral T>
struct Carry_Result
{
    T carry;
    size_t size;
};

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
func add_carry_step(T left, T right, T carry) -> Carried<T>
{
    const T res_low = low_bits(left) + low_bits(right) + carry;
    const T middle_carry = high_bits(res_low);
    
    const T res_high = high_bits(left) + high_bits(right) + middle_carry;
    const T new_carry = high_bits(res_high);

    return Carried<T>{combine_bits(res_low, res_high), new_carry};
}

template <integral T>
func sub_carry_step(T left, T right, T carry) -> Carried<T>
{
    const T res_low = low_bits(left) - low_bits(right) - carry;
    const T middle_carry = res_low >> (BIT_SIZE<T> - 1);

    const T res_high = high_bits(left) - high_bits(right) - middle_carry;
    const T new_carry = res_high >> (BIT_SIZE<T> - 1);

    return Carried<T>{combine_bits(res_low, res_high), new_carry};
}

template <integral T>
func add_carry_patch_step(T left, T carry) -> Carried<T>
{
    T res = left + carry;
    return Carried<T>{res, cast(T)(res < left)};
}

template <integral T>
func sub_carry_patch_step(T left, T carry) -> Carried<T>
{
    T res = left - carry;
    return Carried<T>{res, cast(T)(res > left)};
}

template <integral T>
func add_carry(span<T>* to, span<const T> left, span<const T> right, T base_carry = 0) -> T
{
    assert(to->size() >= right.size());
    assert(to->size() >= left.size());
    if(left.size() < right.size())
        std::swap(left, right);

    T carry = base_carry;
    proc update = [&](size_t i, Carried<T> carried) {
        (*to)[i] = carried.value;
        carry = carried.carry;
    };

    const size_t to_size = min(left.size(), right.size());
    for (size_t i = 0; i < to_size; i++)
        update(i, add_carry_step<T>(left[i], right[i], carry));

    for (size_t i = right.size(); i < left.size(); i++)
        update(i, add_carry_patch_step<T>(left[i], carry));

    return carry;
}

template <integral T>
func sub_carry(span<T>* to, span<const T> left, span<const T> right, T base_carry = 0) -> T
{
    assert(to->size() >= right.size());
    assert(to->size() >= left.size());

    T carry = base_carry;
    proc update = [&](size_t i, Carried<T> carried) {
        (*to)[i] = carried.value;
        carry = carried.carry;
    };

    const size_t to_size = min(left.size(), right.size());
    for (size_t i = 0; i < to_size; i++)
        update(i, sub_carry_step<T>(left[i], right[i], carry));

    for (size_t i = left.size(); i < right.size(); i++)
        update(i, sub_carry_step<T>(0, right[i], carry));

    for (size_t i = right.size(); i < left.size(); i++)
        update(i, sub_carry_patch_step<T>(left[i], carry));

    return carry;
}

template <integral T>
func complement_carry_step(T val, T carry) -> Carried<T>
{
    const T val_inv = ~val;

    const T complement_low = low_bits(val_inv) + carry;
    const T middle_carry = high_bits(complement_low);

    const T complement_high = high_bits(val_inv) + middle_carry;
    const T new_carry = high_bits(complement_high);

    return Carried<T>{combine_bits(complement_low, complement_high), new_carry};
}

template <integral T>
proc twos_complement_carry(span<T>* to, span<const T> left, T base_carry = 1) -> T
{
    assert(to->size() >= left.size());
    assert(is_striped_form(left));

    T carry = base_carry;
    for(size_t i = 0; i < left.size(); i++)
    {
        let res = complement_carry_step<T>(left[i], carry);
        (*to)[i] = res.value;
        carry = res.carry;
    }

    return carry;
}

template <integral T>
func shift_up_carry(span<T>* out, span<const T> in, T by_bits, T prev_high = 0) -> T
{
    assert(out->size() >= in.size());
    size_t shift_items = by_bits / BIT_SIZE<T>;
    size_t shift_bits = by_bits - shift_items * BIT_SIZE<T>;

    for (size_t i = shift_items; i < in.size(); i++)
    {
        T item = in[i - shift_items];
        T low = low_bits(item, shift_bits);
        T high = high_bits(item, shift_bits);

        T composed = prev_high | low << shift_bits;
        
        (*out) = composed;
        prev_high = high;
    }

    null_n(out->data(), min(shift_items, out->size()));
    return prev_high;
}

template <integral T>
func shift_down_carry(span<T>* out, span<const T> in, T by_bits, T prev_prev_lowhigh = 0) -> T
{
    assert(out->size() >= in.size());
    size_t shift_items = by_bits / BIT_SIZE<T>;
    size_t shift_bits = by_bits - shift_items * BIT_SIZE<T>;

    size_t capped_size = min(shift_items, out->size());
    for (size_t i = in.size() - capped_size; i-- > 0;)
    {
        T item = in[i + shift_items];
        T low = low_bits(item, shift_bits);
        T high = high_bits(item, shift_bits);

        T composed = prev_low | high >> shift_bits;

        (*out) = composed;
        prev_low = low;
    }

    null_n(out->data() - capped_size, capped_size);
    return prev_low;
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
func mul_carry_step(T left, T right, T carry) -> Carried<T>
{
    //we do the carried multiplication by multiplying each digit normally and summing the overflow 
    // => 25  
    //    *5
    //  ----
    //    25 <- here 2 is 'carried over' to the next iteration 
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
    //  => s_mixed slot is half in the next iteration => we dont shift (cause we would shift out of type)
    //  => s_next slot is entirely in the next iteration => also no shift
    //
    // s_mixed = s_mixed_ns * b   //ns means 'no shift'
    // s_next = s_next_ns * b 
    size_t h = HALF_BIT_SIZE<T>;

    T l = left;
    T r = right;

    T l1 = low_bits(l);
    T l2 = high_bits(l);
    T r1 = low_bits(r);
    T r2 = high_bits(r);

    T s_cur = l1*r1;
    T s_mixed_ns = l1*r2 + l2*r1; 
    T s_next_ns = l2*r2;

    T curr_value = s_cur + (low_bits(s_mixed_ns) << h); // even though we are taking lower half of s_mixed we have to insert it to upper half
    T next_value = s_next_ns + high_bits(s_mixed_ns);

    const T carried = curr_value + last_value;

    //if overflowed next carry is 1 bigger
    if(carried < curr_value)
         next_value += 1;

    return Carried<T>{carried, next_value};
}

template <integral T>
func is_power_of_two(T num)
{
    return num != 0 && (num & (num - 1)) == 0;
}

static constexpr bool DO_MUL_OPTIMS = true;
template <integral T>
func mul_carry(span<T>* to, span<const T> left, T right) -> Carry_Result<T>
{
    //assert(is_striped_form(left));
    assert(to->size() >= left.size());

    if(right == 0) 
        return Carry_Result<T>{0, 0};

    if constexpr(DO_MUL_OPTIMS)
    {
        if(right == 1)
        {
            copy_n(to->data(), left.data(), left.size());
            return Carry_Result<T>{0, left.size()};
        }

        if(is_power_of_two(right))
        {
            let res = shift_up_carry(to, left, right);
            return Carry_Result<T>{res, left.size()};
        }
    }

    T last_value = 0;
    for (size_t i = 0; i < left.size(); i++)
    {
        let res = mul_carry_step<T>(left[i], right, last_value);
        (*to)[i] = res.value;
        last_value = res.carry;
    }

    let obtained = as_number(*to);

    return Carry_Result<T>{
        last_value,
        left.size(),
    };
}

template <integral T>
func mul_quadratic(span<T>* to, span<T>* temp, span<const T> left, span<const T> right)
{
    size_t max_result_size = left.size() + right.size() + 1;
    assert(to->size() >= max_result_size);
    assert(temp->size() >= right.size() + 1);
    //assert(is_striped_form(left));
    //assert(is_striped_form(right));
    null_n(to->data(), to->size());

    for(size_t i = 0; i < right.size(); i++)
    {
        let mul_res = mul_carry<T>(temp, left, right[i]);
        size_t curr_size = mul_res.size + 1;
        (*temp)[mul_res.size] = mul_res.carry;
        
        mut shifted_to = to->subspan(i);
        assert(shifted_to.size() >= curr_size);
        assert(temp->size() >= curr_size);

        let carry = add_carry<T>(&shifted_to, shifted_to, temp->subspan(0, curr_size));
        assert(carry == 0); //no carry should happen in this case
    }
}

template <integral T>
func mul_quadratic_fused(span<T>* to, span<const T> left, span<const T> right)
{
    size_t max_result_size = left.size() + right.size() + 1;
    assert(to->size() >= max_result_size);
    null_n(to->data(), to->size());

    for(size_t i = 0; i < right.size(); i++)
    {
        let mul_res = mul_carry<T>(temp, left, right[i]);
        size_t curr_size = mul_res.size + 1;
        (*temp)[mul_res.size] = mul_res.carry;

        assert(shifted_to.size() >= curr_size);
        assert(temp->size() >= curr_size);

        let carry = add_carry<T>(&shifted_to, shifted_to, temp->subspan(0, curr_size));
        assert(carry == 0); //no carry should happen in this case

        //MUL
        mut shifted_to = to->subspan(i);
        T curr_right = right[i];
        if(curr_right == 0) 
            continue;

        if constexpr(DO_MUL_OPTIMS)
        {
            if(curr_right == 1)
            {
                let carry = add_carry(&shifted_to, shifted_to, left);
                assert(carry == 0); //no carry should happen in this case
                continue;
            }

            //skip for now
            /*if(is_power_of_two(curr_right))
            {
                let res = shift_up_carry(to, left, curr_right);
                return Carry_Result<T>{res, left.size()};
            }*/
        }

        T last_mul_val = 0;
        for (size_t i = 0; i < left.size(); i++)
        {
            T curr_left = left[i];
            let res = mul_carry_step<T>(curr_left, curr_right, last_mul_val);


            last_mul_val = res.carry;
        }

        let obtained = as_number(*to);

        return Carry_Result<T>{
            last_value,
                left.size(),
        };

        //ADD
        T carry = base_carry;
        proc update = [&](size_t i, Carried<T> carried) {
            (*to)[i] = carried.value;
            carry = carried.carry;
        };

        const size_t to_size = min(left.size(), right.size());
        for (size_t i = 0; i < to_size; i++)
            update(i, add_carry_step<T>(left[i], right[i], carry));

        for (size_t i = right.size(); i < left.size(); i++)
            update(i, add_carry_patch_step<T>(left[i], carry));

    }
}
