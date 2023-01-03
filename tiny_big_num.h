#pragma once

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <type_traits>

namespace tiny_num
{
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

    #define let const auto
    #define mut auto
    #define proc constexpr auto
    #define func [[nodiscard]] constexpr auto
    #define cast(...) (__VA_ARGS__)

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
    func slice(Slice<T> items, size_t from) -> Slice<T> {
        return slice_size<T>(items, from, items.size - from);
    }

    template <typename T>
    func trim(Slice<T> items, size_t to_size) -> Slice<T> {
        return slice_size<T>(items, 0, to_size);
    }

    template<typename T>
    struct Id {
        using type = T;
    };

    template<typename T>
    using No_Infer = Id<T>::type;

    template<typename T>
    func max(T a, No_Infer<T> b) -> T
    {
        return a > b ? a : b;
    }

    template<typename T>
    func min(T a, No_Infer<T> b) -> T
    {
        return a < b ? a : b;
    }

    template<typename T>
    func div_round_up(T value, No_Infer<T> to_multiple_of) -> T
    {
        static_assert(std::is_integral_v<T>, "!");
        return (value + to_multiple_of - 1) / to_multiple_of;
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

    //template<typename T>
    //static constexpr bool is_allowed_digit = std::is_unsigned_v<T>;

    template<typename T>
    static constexpr bool is_allowed_digit = true;

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
        static_assert(is_allowed_digit<T>, "!");
        return cast(T) (FULL_MASK<T> << index);
    }

    template <typename T>
    func low_mask(size_t index = HALF_BIT_SIZE<T>) -> T {
        static_assert(is_allowed_digit<T>, "!");
        return cast(T) ~high_mask<T>(index);
    }

    template <typename T>
    func high_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
        static_assert(is_allowed_digit<T>, "!");
        assert(index < BIT_SIZE<T>);
        return value >> index;
    }

    template <typename T>
    func low_bits(T value, size_t index = HALF_BIT_SIZE<T>) -> T {
        static_assert(is_allowed_digit<T>, "!");
        let mask = low_mask<T>(index);
        return value & mask;
    };  

    template <typename T>
    func combine_bits(T low, T high, size_t index = HALF_BIT_SIZE<T>) -> T {
        assert(index < BIT_SIZE<T>);
        static_assert(is_allowed_digit<T>, "!");
        return low_bits(low, index) | (high << index);
    };

    template <typename T>
    func dirty_combine_bits(T low, T high, size_t index = HALF_BIT_SIZE<T>) -> T {
        assert(index < BIT_SIZE<T>);
        static_assert(is_allowed_digit<T>, "!");
        assert(high_bits(low, index) == 0 && "low must not have high bits use combine_bits instead");
        return low | (high << index);
    };

    template <typename T>
    func find_last_set_bit(T val) -> T  
    {
        static_assert(is_allowed_digit<T>, "!");
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

        static_assert(is_allowed_digit<T>, "!");
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
        static_assert(is_allowed_digit<T>, "!");
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
            size_t count = 0;
            for (; val != 0; val >>= 1)
                if (val & 1)
                    count++;

            return count;
        }
    }

    template <typename T>
    struct Single_Overflow
    {
        T value = 0;
        T overflow = 0;

        bool constexpr operator ==(Single_Overflow const&) const noexcept = default;
    };

    template <typename T>
    struct Batch_Overflow
    {
        Slice<T> slice;
        T overflow = 0;

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
        static_assert(is_allowed_digit<T>, "!");

        size_t i = num.size;
        for(; i-- > 0;)
            if(num[i] != 0)
                break;

        return i;
    }

    template <typename T>
    func find_first_set_digit(Slice<T> num) -> size_t
    {
        static_assert(is_allowed_digit<T>, "!");
        
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

    template <typename Digit, typename Number>
    func digits_to_represent() -> size_t
    {
        return (sizeof(Number) + sizeof(Digit) - 1) / sizeof(Digit);
    }

    template <typename To = std::uint64_t, typename T = size_t>
    func to_number(Slice<T> bignum) -> To
    {
        static_assert(is_allowed_digit<T>, "!");
        let stripped = striped_trailing_zeros(bignum);

        const size_t total_size = stripped.size * sizeof(T);
        assert(total_size <= sizeof(To));

        To out = 0;
        for(size_t i = 0; i < stripped.size; i++)
        {
            size_t shift_by = i * BIT_SIZE<T>;
            assert(shift_by < BIT_SIZE<To>);
            out |= cast(To) stripped[i] << shift_by;
        }

        return out;
    }

    template <typename T, typename From = std::uint64_t>
    proc from_number(Slice<T>* bignum, From from) -> Slice<T>
    {
        static_assert(is_allowed_digit<T>, "!");
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
    namespace detail
    {
        template <typename T, typename... Ts>
        func single_add_overflow_any(T carry, Ts... args) -> Single_Overflow<T>
        {
            static_assert(is_allowed_digit<T>, "!");
            static_assert((std::is_same_v<T, Ts> && ...), "!");
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
            constexpr size_t arg_count = sizeof...(args) + 1;
            constexpr size_t max_args = cast(size_t) (cast(T) 1 << HALF_BIT_SIZE<T>);
            static_assert(arg_count <= max_args && "too many operands");

            const T res_low = (low_bits(cast(T) args) + ...) + carry;
            const T middle_carry = high_bits(res_low);

            const T res_high = (high_bits(cast(T) args) + ...) + middle_carry;
            const T new_carry = high_bits(res_high);

            return Single_Overflow<T>{combine_bits(res_low, res_high), new_carry};
        }

        template <typename T, typename... Ts>
        func single_sub_overflow_any(T carry, T left, Ts... args) -> Single_Overflow<T>
        {
            static_assert(is_allowed_digit<T>, "!");
            static_assert((std::is_same_v<T, Ts> && ...), "!");

            constexpr size_t arg_count = sizeof...(args) + 1;
            constexpr size_t max_args = cast(size_t) (cast(T) 1 << HALF_BIT_SIZE<T>);
            static_assert(arg_count < max_args && "too many operands");

            const T res_low = high_mask<T>() + low_bits(left) - (low_bits(args) + ...) - carry;
            const T middle_carry = low_mask<T>() - high_bits<T>(res_low);

            const T res_high = high_mask<T>() + high_bits(left) - (high_bits(args) + ...) - middle_carry;
            const T new_carry = low_mask<T>() - high_bits<T>(res_high);

            return Single_Overflow<T>{combine_bits(res_low, res_high), new_carry};
        }

    }

    template <typename T, typename... Ts>
    func single_add_no_overflow(T first, Ts... args) -> T
    {
        static_assert((std::is_same_v<T, Ts> && ...), "!");
        #ifdef NDEBUG
            return first + (args + ...);
        #else
            const size_t size = 1 + sizeof...(args);
            const T arr[] = {first, args...};
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

    template <typename T>
    func single_add_overflow(T left, T right, T carry = 0) -> Single_Overflow<T> {
        return detail::single_add_overflow_any<T>(carry, left, right);
    }

    template <typename T>
    func single_sub_overflow(T left, T right, T carry = 0) -> Single_Overflow<T> {
        return detail::single_sub_overflow_any<T>(carry, left, right);
    }

    template <typename T, typename... Ts>
    func single_add_overflow_any(T left, T right, Ts... rest) -> Single_Overflow<T> {
        return detail::single_add_overflow_any<T>(0, left, right, rest...);
    }

    template <typename T, typename... Ts>
    func single_sub_overflow_any(T left, T right, Ts... rest) -> Single_Overflow<T> {
        return detail::single_sub_overflow_any<T>(0, left, right, rest...);
    }

    enum class Location
    {
        IN_PLACE,
        OUT_OF_PLACE
    };

    template <typename T>
    func batch_add_or_sub_overflow_short(Slice<T>* to, Slice<const T> left, T carry, bool is_addition, Location location = Location::OUT_OF_PLACE, size_t from = 0) -> Batch_Overflow<T>
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
            copy_n<T>(to->data + j, left.data + j, left.size - j, Iter_Direction::FORWARD);

        return Batch_Overflow<T>{trimmed_to, carry};
    }

    //add short overflow
    template <typename T>
    func batch_add_overflow_short(Slice<T>* to, Slice<const T> left, T right, Location location = Location::OUT_OF_PLACE, size_t from = 0) -> Batch_Overflow<T>
    {
        return batch_add_or_sub_overflow_short(to, left, right, true, location, from);
    }

    template <typename T>
    func batch_sub_overflow_short(Slice<T>* to, Slice<const T> left, T right, Location location = Location::OUT_OF_PLACE, size_t from = 0) -> Batch_Overflow<T>
    {
        return batch_add_or_sub_overflow_short(to, left, right, false, location, from);
    }

    //add long overflow
    template <typename T>
    func batch_add_overflow_long(Slice<T>* to, Slice<const T> left, Slice<const T> right, Location location = Location::OUT_OF_PLACE, T carry_in = 0) -> Batch_Overflow<T>
    {
        static_assert(is_allowed_digit<T>, "!");
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
            let res = single_add_overflow<T>(left[i], right[i], carry);
            (*to)[i] = res.value;
            carry = res.overflow;
        }

        return batch_add_overflow_short<T>(to, left, carry, location, right.size);
    }

    template <typename T>
    func batch_sub_overflow_long(Slice<T>* to, Slice<const T> left, Slice<const T> right, Location location = Location::OUT_OF_PLACE, T carry_in = 0) -> Batch_Overflow<T>
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(to->size >= right.size);
        assert(to->size >= left.size);
        assert(check_are_one_way_aliasing<T>(left, *to) == false);
        assert(carry_in == 0 || carry_in == 1);
        //assert(high_bits(carry_in) == 0);

        T carry = carry_in;
        proc update = [&](size_t i, Single_Overflow<T> oveflow) {
            (*to)[i] = oveflow.value;
            carry = oveflow.overflow;
        };

        const size_t min_size = min(left.size, right.size);
        const size_t max_size = max(left.size, right.size);
        for (size_t i = 0; i < min_size; i++)
            update(i, single_sub_overflow<T>(left[i], right[i], carry));

        for (size_t i = left.size; i < right.size; i++)
            update(i, single_sub_overflow<T>(0, right[i], carry));

        if(right.size < left.size)
            return batch_sub_overflow_short<T>(to, left, carry, location, right.size);

        return Batch_Overflow<T>{trim(*to, max_size), carry};
    }

    template <typename T>
    func single_complement_overflow(T val, T carry) -> Single_Overflow<T>
    {
        const T val_inv = ~val;

        const T complement_low = low_bits(val_inv) + carry;
        const T middle_carry = high_bits(complement_low);

        const T complement_high = high_bits(val_inv) + middle_carry;
        const T new_carry = high_bits(complement_high);

        return Single_Overflow<T>{combine_bits(complement_low, complement_high), new_carry};
    }

    template <typename T>
    proc batch_complement_overflow(Slice<T>* to, Slice<const T> left, T carry_in = 1) -> Batch_Overflow<T>
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(to->size >= left.size);
        assert(is_striped_number(left));
        assert(check_are_one_way_aliasing<T>(left, *to) == false);
        assert(carry_in == 1 || carry_in == 0);

        T carry = carry_in;
        for(size_t i = 0; i < left.size; i++)
        {
            let res = single_complement_overflow<T>(left[i], carry);
            (*to)[i] = res.value;
            carry = res.overflow;
        }

        return Batch_Overflow<T>{trim(*to, left.size), carry};
    }

    template <typename T>
    func single_shift_up_overflow(T low_item, size_t by_bits, T high_item) -> Single_Overflow<T>
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
        return Single_Overflow<T>{composed, shifted_out};
    }

    template <typename T>
    func single_shift_down_overflow(T low_item, size_t by_bits, T high_item) -> Single_Overflow<T>
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
        return Single_Overflow<T>{composed, shifted_out};
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
    proc batch_shift_up_overflow(Slice<T>* out, Slice<const T> in, size_t by_bits, 
        Iter_Direction direction = Iter_Direction::FORWARD, T carry_in = 0) -> Batch_Overflow<T>
    {
        static_assert(is_allowed_digit<T>, "!");
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
                let res = single_shift_up_overflow<T>(prev, by_bits, curr);
                (*out)[i] = res.value;
                prev = curr;
            }

            let carry_out = single_shift_up_overflow<T>(prev, by_bits, 0);
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
                (*out)[i] = single_shift_up_overflow<T>(curr, by_bits, prev).value;
                prev = curr;
            }

            {
                let res = single_shift_up_overflow<T>(carry_in, by_bits, prev);
                (*out)[0] = res.value;
            }


            let carry_out = single_shift_up_overflow<T>(shifted_out, by_bits, 0);
            return Batch_Overflow<T>{out_trimmed, carry_out.value};
        }
    }

    template <typename T>
    proc batch_shift_down_overflow(Slice<T>* out, Slice<const T> in, size_t by_bits, 
        Iter_Direction direction = Iter_Direction::FORWARD, T carry_in = 0) -> Batch_Overflow<T>
    {
        static_assert(is_allowed_digit<T>, "!");
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
                let res = single_shift_down_overflow<T>(prev, by_bits, curr);
                (*out)[i] = res.value;
                prev = curr;
            }

            {
                let res = single_shift_down_overflow<T>(prev, by_bits, carry_in);
                (*out)[in.size - 1] = res.value;
            }

            let carry_out = single_shift_down_overflow<T>(0, by_bits, shifted_out);
            return Batch_Overflow<T>{out_trimmed, carry_out.value};
        }
        else
        {
            assert(check_are_one_way_aliasing<T>(*out, in) == false);

            T prev = carry_in;
            for (size_t i = in.size; i-- > 0;)
            {
                let curr = in[i];
                (*out)[i] = single_shift_down_overflow<T>(curr, by_bits, prev).value;
                prev = curr;
            }

            let carry_out = single_shift_down_overflow<T>(0, by_bits, prev);
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
    func single_mul_overflow(T left, T right, T last_value, Mul_Overflow_Optims const& mul_optims = Mul_Overflow_Optims::NONE) -> Single_Overflow<T>
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
            let s_mixed_ns_res = single_add_overflow<T>(l1*r2, l2*r1);

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

        let curr_value = single_add_overflow_any<T>(s_cur, low_mixed, last_value); //also add last_value to save ops
        T next_value =   single_add_no_overflow<T>(s_next_ns, high_mixed, curr_value.overflow); //also add the overflow

        return Single_Overflow<T>{curr_value.value, next_value};
    }

    template <typename T>
    func single_is_power_of_two(T num)
    {
        return num != 0 && (num & (num - 1)) == 0;
    }

    template <typename T>
    proc batch_mul_overflow(Slice<T>* to, Slice<const T> left, T right, Optims const& optims, T carry_in = 0) -> Batch_Overflow<T>
    {
        static_assert(is_allowed_digit<T>, "!");
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
            if(single_is_power_of_two(right))
            {
                //@NOTE: here carry is handled differently the in the rest of the mul algorhitm
                //  this doesnt affect how the algorhitm appears from the outside not does it affect the
                //  ability for batches to be chained together with carries but it forbids any kind of dependence
                //  on the carry in/out in this case. This is not ideal and should be changed so that the
                //  shift itself handles carry in similar way to mul
                let shift_bits = find_last_set_bit(right);
                return batch_shift_up_overflow<T>(to, left, shift_bits, Iter_Direction::FORWARD, carry_in);
            }
        }

        #define do_carry_loop(operation)            \
        {                                           \
            for (size_t i = 0; i < left.size; i++)  \
            {                                       \
                const Single_Overflow<T> res = operation;  \
                (*to)[i] = res.value;               \
                carry = res.overflow;               \
            }                                       \
        }                                           \

        if(optims.mul_half_bits && high_bits(right) == 0)
            do_carry_loop((single_mul_overflow<T>(left[i], right, carry, Mul_Overflow_Optims::LOW_BITS_ONLY)))
        else if(optims.mul_half_bits && low_bits(right) == 0)
            do_carry_loop((single_mul_overflow<T>(left[i], right, carry, Mul_Overflow_Optims::HIGH_BITS_ONLY)))
        else
            do_carry_loop((single_mul_overflow<T>(left[i], right, carry)))

        #undef do_carry_loop
        return Batch_Overflow<T>{trimmed_to, carry};
    }

    template <typename T>
    func single_div_overflow(T left, T right, T carry_in = 0) -> Single_Overflow<T>
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

        const T operand_high = combine_bits<T>(high_bits(left), carry_in);
        const T res_high = operand_high / right;
        const T middle_carry = operand_high % right;

        const T operand_low = combine_bits<T>(low_bits(left), middle_carry);
        const T res_low = operand_low / right;
        const T out_carry = operand_low % right;

        const T res = dirty_combine_bits(res_low, res_high);
        return Single_Overflow<T>{res, out_carry};
    }

    template <typename T>
    proc batch_div_overflow(Slice<T>* to, Slice<const T> left, T right, Optims const& optims, T carry_in = 0) -> Batch_Overflow<T>
    {
        static_assert(is_allowed_digit<T>, "!");
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
            if(single_is_power_of_two(right))
            {
                let shift_bits = find_last_set_bit(right);
                let shift_res = batch_shift_down_overflow<T>(to, left, shift_bits, Iter_Direction::BACKWARD, carry_in);
                //because of the way shifting down result is defined we have to process the reuslt
                const T carry = shift_res.overflow >> (BIT_SIZE<T> - shift_bits);
                return Batch_Overflow<T>{trimmed_to, carry};
            }
        }

        T carry = carry_in;
        for (size_t i = left.size; i-- > 0; )
        {
            let res = single_div_overflow<T>(left[i], right, carry);
            (*to)[i] = res.value;
            carry = res.overflow;
        }

        return Batch_Overflow<T>{trimmed_to, carry};
    }

    template <typename T>
    proc batch_rem_overflow(Slice<const T> left, T right, Optims const& optims, T carry_in = 0) -> T
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(high_bits(right) == 0 && "only works for divisors under half bit size");
        assert(right != 0 && "cannot divide by zero");

        if(optims.rem_optims)
        {
            if(single_is_power_of_two(right))
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
            let res = single_div_overflow<T>(left[i], right, carry);
            carry = res.overflow;
        }

        return carry;
    }


    func required_add_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return max(left_size, right_size) + 1;
    }

    template <typename T>
    proc add(Slice<T>* to, Slice<const T> left, Slice<const T> right, Location location = Location::OUT_OF_PLACE)
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(is_striped_number(left));
        assert(is_striped_number(right));
        assert(to->size >= required_add_out_size(left.size, right.size));

        let res = batch_add_overflow_long<T>(to, left, right, location);
        if(res.overflow == 0)
        {
            if(location == Location::IN_PLACE)
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

    func required_sub_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return left_size;
    }

    template <typename T>
    proc sub(Slice<T>* to, Slice<const T> left, Slice<const T> right, Location location = Location::OUT_OF_PLACE)
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(is_striped_number(left));
        assert(is_striped_number(right));
        assert(to->size >= required_sub_out_size(left.size, right.size));

        let res = batch_sub_overflow_long<T>(to, left, right, location);
        assert(res.overflow == 0 && "left must be bigger than right");

        return striped_trailing_zeros(res.slice);
    }

    func required_mul_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return left_size + right_size;
    }

    template <typename T>
    proc mul(Slice<T>* to, Slice<const T> left, T right, Optims const& optims)
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(is_striped_number(left));
        assert(to->size >= required_sub_out_size(left.size, 1));

        let res = batch_mul_overflow<T>(to, left, right, optims);
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
    struct Div_Result
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
    proc div_bit_by_bit(Slice<T>* quotient, Slice<T>* remainder, Slice<const T> num, Slice<const T> den, Optims const& optims) -> Div_Result<T>
    {
        static_assert(is_allowed_digit<T>, "!");
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

        assert(den.size > 0 && "cannot divide by zero!");

        if(num.size == 0)
        {
            let stripped_quotient = trim(*quotient, 0);
            let stripped_remainder = trim(*remainder, 0);
            return Div_Result<T>{stripped_quotient, stripped_remainder};
        }

        if(num.size < den.size)
        {
            assert(remainder->size >= num.size);
            mut stripped_remainder = trim(*remainder, num.size);
            mut stripped_quotient = trim(*quotient, 0);

            copy_slice<T>(&stripped_remainder, num, Iter_Direction::NO_ALIAS);
            return Div_Result<T>{stripped_quotient, stripped_remainder};
        }

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
                let div_res = batch_div_overflow<T>(&trimmed_quotient, num, last(den), optims);
                stripped_quotient = striped_trailing_zeros<T>(div_res.slice);
                remainder_val = div_res.overflow;
            }
            else
            {
                remainder_val = batch_rem_overflow<T>(num, last(den), optims);
                stripped_quotient = trim(trimmed_quotient, 0);
            }

            size_t remainder_size = 0;
            if(remainder_val != 0)
            {
                trimmed_remainder[0] = remainder_val;
                remainder_size = 1;
            }

            let stripped_remainder = trim(trimmed_remainder, remainder_size);
            return Div_Result<T>{stripped_quotient, stripped_remainder};
        }

        Slice<T> curr_remainder = trim(trimmed_remainder, 0);
        for(size_t i = BIT_SIZE<T> * num.size; i-- > 0; ) 
        {
            let last_size = curr_remainder.size;
            let shift_res = batch_shift_up_overflow<T>(&curr_remainder, curr_remainder, 1);
            const T num_ith_bit = get_nth_bit(num, i);

            if(shift_res.overflow != 0 || (last_size == 0 && num_ith_bit == 1))
            {
                curr_remainder = trim(trimmed_remainder, last_size + 1);
                curr_remainder[last_size] = shift_res.overflow;
            }

            set_nth_bit(&trimmed_remainder, 0, num_ith_bit);
            if(compare<T>(curr_remainder, den) >= 0)
            {
                let sub_res = batch_sub_overflow_long<T>(&curr_remainder, curr_remainder, den, Location::IN_PLACE);
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

        return Div_Result<T>{stripped_quotient, curr_remainder};
    }

    template <typename T>
    proc div(Slice<T>* quotient, Slice<T>* remainder, Slice<const T> num, Slice<const T> den, Optims const& optims) -> Div_Result<T>
    {
        return div_bit_by_bit(quotient, remainder, num, den, optims);
    }

    template <typename T>
    proc rem(Slice<T>* remainder, Slice<const T> num, Slice<const T> den, Optims const& optims) -> Slice<T>
    {
        static_assert(is_allowed_digit<T>, "!");
        Slice<T> quotient = {nullptr, 0};
        let div_res = div_bit_by_bit<T, false>(&quotient, remainder, num, den, optims);

        assert(div_res.quotient.size == 0 && "in case of only remainder the quotient size should be 0");
        assert(div_res.quotient.data == nullptr && "and data nullptr as we set it");
        return div_res.remainder;
    }

    template <typename T>
    proc batch_fused_mul_add_overflow(Slice<T>* to, Slice<const T> added, Slice<const T> multiplied, T coeficient, 
        Optims const& optims, const Location location = Location::OUT_OF_PLACE, T add_carry = 0, T mul_carry = 0) -> Batch_Overflow<T>
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
                if (location == Location::IN_PLACE)
                    return Batch_Overflow<T>{trim(*to, 0), 0};
                else
                {
                    copy_n(to->data, added.data, added.size, Iter_Direction::FORWARD);
                    return Batch_Overflow<T>{trimmed_to, 0};
                }
            }

            if(coeficient == 1)
                return batch_add_overflow_long<T>(to, added, multiplied, location);
        }

        Slice<T> trimmed_to = trim(*to, added.size);
        if(optims.mul_shift 
            && optims.mul_consts //zero check must be performed (else shifting is UB!)
            && single_is_power_of_two(coeficient))
        {
            let shift_bits = find_last_set_bit(coeficient);
            const Slice<T> trimmed_to = trim(*to, added.size);
            T prev_muled = mul_carry;
            for (size_t j = 0; j < multiplied.size; j++)
            {
                const T curr_muled = multiplied[j];
                const T curr_added = added[j];

                let mul_res = single_shift_up_overflow<T>(prev_muled, shift_bits, curr_muled);
                let add_res = single_add_overflow<T>(curr_added, mul_res.value, add_carry);

                trimmed_to[j] = add_res.value;

                prev_muled = curr_muled;
                add_carry = add_res.overflow;
            }

            mul_carry = single_shift_up_overflow<T>(prev_muled, shift_bits, 0).value;
        }
        else
        {
            for (size_t j = 0; j < multiplied.size; j++)
            {
                T curr_muled = multiplied[j];
                T curr_added = added[j];

                let mul_res = single_mul_overflow<T>(curr_muled, coeficient, mul_carry);
                let add_res = single_add_overflow<T>(curr_added, mul_res.value, add_carry);

                trimmed_to[j] = add_res.value;

                mul_carry = mul_res.overflow;
                add_carry = add_res.overflow;
            }
        }
        
        const T combined_carry = single_add_no_overflow(mul_carry, add_carry);
        return batch_add_overflow_short<T>(&trimmed_to, added, combined_carry, location, multiplied.size);
    }

    template <typename T>
    proc mul(Slice<T>* to, Slice<T>* aux, Slice<const T> left, Slice<const T> right, Optims const& optims, size_t depth = 0);

    func required_mul_quadratic_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return required_mul_out_size(left_size, right_size);
    }

    func required_mul_quadratic_aux_size(size_t left_size, size_t right_size) -> size_t
    {
        return required_mul_out_size(left_size, right_size);
    }

    template <typename T>
    proc mul_quadratic(Slice<T>* to, Slice<T>* temp, Slice<const T> left, Slice<const T> right, Optims const& optims) -> Slice<T>
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(check_are_aliasing<T>(*to, left) == false);
        assert(check_are_aliasing<T>(*to, right) == false);
        assert(check_are_aliasing<T>(*temp, left) == false);
        assert(check_are_aliasing<T>(*temp, right) == false);
        assert(check_are_aliasing<T>(*to, *temp) == false);

        if(left.size < right.size)
            std::swap(left, right);

        const size_t max_result_size = required_mul_out_size(left.size, right.size);
        const size_t max_temp_size = required_mul_out_size(left.size, right.size);
        assert(to->size >= max_result_size);
        assert(temp->size >= max_temp_size);
        
        Slice<T> trimmed_to = trim(*to, max_result_size);
        Slice<T> trimmed_temp = trim(*temp, max_temp_size);
        null_n(trimmed_to.data, trimmed_to.size);

        for(size_t i = 0; i < right.size; i++)
        {
            let mul_res = batch_mul_overflow<T>(&trimmed_temp, left, right[i], optims);

            //batch_mul_overflow is limited to left.size elems => in case of overflow
            // even if the overflow would fit into temp it is not added => we add it back
            let mul_size = mul_res.slice.size;
            let mul_slice = trim(trimmed_temp, mul_size + 1); 
            mul_slice[mul_size] = mul_res.overflow;

            mut shifted_to = slice(trimmed_to, i);
            assert(shifted_to.size >= mul_slice.size);

            let add_res = batch_add_overflow_long<T>(&shifted_to, shifted_to, mul_slice, Location::IN_PLACE);
            assert(add_res.overflow == 0 && "no final carry should be left");
        }

        return striped_trailing_zeros(trimmed_to);
    }


    func required_mul_quadratic_fused_out_size(size_t left_size, size_t right_size) -> size_t
    {
        return required_mul_out_size(left_size, right_size);
    }

    //@TODO: make it so that to is only operated on to the required sized
    template <typename T>
    proc mul_quadratic_fused(Slice<T>* to, Slice<const T> left, Slice<const T> right, Optims const& optims) -> Slice<T>
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(check_are_aliasing<T>(*to, left) == false);
        assert(check_are_aliasing<T>(*to, right) == false);

        if(left.size < right.size)
            std::swap(left, right);

        const size_t max_result_size = required_mul_out_size(left.size, right.size);
        assert(to->size >= max_result_size);
        Slice<T> trimmed_to = trim(*to, max_result_size);
        null_n(trimmed_to.data, trimmed_to.size);

        //Slice<T> last_modified = slice(trimmed_to, 0); //@TODO
        for(size_t i = 0; i < right.size; i++)
        {
            mut shifted_to = slice(trimmed_to, i);

            let mul_add_res = batch_fused_mul_add_overflow<T>(&shifted_to, shifted_to, left, right[i], optims, Location::IN_PLACE);

            assert(mul_add_res.overflow == 0 && "no carry should happen in this case");
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
            std::swap(left_size, right_size);

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


    template <typename T>
    proc mul_karatsuba(Slice<T>* to, Slice<T>* aux, Slice<const T> left, Slice<const T> right, Optims const& optims, size_t depth = 0, bool is_run_alone = true) -> Slice<T>
    {
        static_assert(is_allowed_digit<T>, "!");

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


        size_t required_out_size = required_mul_out_size(x.size, y.size);
        Slice<T> trimmed_out = trim(*to, required_out_size);

        size_t required_add_x_sum_size = required_add_out_size(x1.size, x2.size);
        size_t required_add_y_sum_size = required_add_out_size(y1.size, y2.size);

        size_t required_z1_size = required_mul_out_size(x1.size, y1.size);
        size_t required_z2_size = required_mul_out_size(required_add_x_sum_size, required_add_y_sum_size); //after multiplies will be even less
        size_t required_z3_size = required_mul_out_size(x2.size, y2.size);

        size_t total_used_aux_size = required_add_x_sum_size + required_add_y_sum_size + required_z2_size;
        size_t returned_required_aux = required_mul_karatsuba_aux_size(left.size, right.size, 0);
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
        Slice<T> x_sum_y_sum_m_z1 = sub<T>(&x_sum_y_sum, x_sum_y_sum, z1, Location::IN_PLACE);
        Slice<T> z2 = sub<T>(&x_sum_y_sum_m_z1, x_sum_y_sum_m_z1, z3, Location::IN_PLACE);

        //instead of multiplying z2 by base we add it to the appropriate position
        Slice<T> out_z2_up = slice(trimmed_out, base);

        let add_res = batch_add_overflow_long<T>(&out_z2_up, out_z2_up, z2, Location::IN_PLACE);
        assert(add_res.overflow == 0 && "should not overflow");

        Slice<T> out = striped_trailing_zeros(trimmed_out);

        return out;
    }

    template <typename T>
    proc mul(Slice<T>* to, Slice<T>* aux, Slice<const T> left, Slice<const T> right, Optims const& optims, size_t depth)
    {
        static_assert(is_allowed_digit<T>, "!");
        assert(is_striped_number(left));
        assert(is_striped_number(right));
        assert(are_aliasing<T>(*to, left) == false);
        assert(are_aliasing<T>(*to, right) == false);
        assert(to->size >= required_mul_out_size(left.size, right.size));

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

    template <typename T>
    func log2(Slice<T> left) -> size_t
    {
        assert(is_striped_number(left));
        if(left.size == 0)
            return 0;

        size_t last_log = find_last_set_bit(last(left));
        return last_log + (left.size - 1) * BIT_SIZE<T>;
    }

    func required_pow_out_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
    {
        if(power == 0 || num_bit_size == 0)
            return 1;

        size_t bit_size = (num_bit_size + 1) * power;
        size_t item_size = div_round_up(bit_size, single_digit_bit_size) + 1;

        return item_size;
    }

    func required_pow_by_squaring_single_aux_swap_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
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

    func required_pow_by_squaring_aux_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
    {
        if(power == 0)
            return 0;
        
        size_t square = required_pow_by_squaring_single_aux_swap_size(num_bit_size, power, single_digit_bit_size);
        size_t out_swap = required_pow_out_size(num_bit_size, power, single_digit_bit_size);

        return square * 2 + out_swap;
    }

    func optimal_pow_aux_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) {
        return required_pow_by_squaring_aux_size(num_bit_size, power, single_digit_bit_size);
    }

    template <typename T>
    func single_pow_trivial(T x, T n) -> T
    {
        T res = 1;
        for(T i = 0; i < n; i++)
            res *= x;

        return res;
    }

    template <typename T>
    func single_pow_by_squaring(T x, T n) -> T
    {
        if (x <= 1) 
            return x;
        if(n == 0)  
            return 1;

        T i = 0;
        T y = 1;
        size_t max_pos = find_last_set_bit(n);
        for(; i <= max_pos; i++)
        {
            T bit = get_bit<T>(n, i);
            if(bit)
                y *= x;
            x *= x;
        }
        return y;
    }

    template <typename T>
    proc pow_by_squaring(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t power, Optims const& optims) -> Slice<T>
    {
        assert(is_striped_number(num));
        //no two can alias
        assert(are_aliasing<T>(*to, num) == false);
        assert(are_aliasing<T>(*aux, num) == false);
        assert(are_aliasing<T>(*aux, *to) == false);

        size_t bit_size = log2(num);
        const size_t required_size = required_pow_out_size(bit_size, power, BIT_SIZE<T>);
        const size_t required_sigle_aux = required_pow_by_squaring_single_aux_swap_size(bit_size, power, BIT_SIZE<T>);

        assert(to->size >= required_size);
        assert(aux->size >= required_pow_by_squaring_aux_size(bit_size, power, BIT_SIZE<T>));

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

    func required_pow_trivial_aux_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
    {
        if(power <= 1)
            return 0;

        //minus one since in the algorhitm we always multiply from output buffer to auxiliary
        // so that the final mutliply will be to the output buffer
        // => the maximum power that will be stored in the auxiliary is one less
        return required_pow_out_size(num_bit_size, power - 1, single_digit_bit_size);
    }
    
    func required_pow_aux_size(size_t num_bit_size, size_t power, size_t single_digit_bit_size) -> size_t
    {
        return required_pow_trivial_aux_size(num_bit_size, power, single_digit_bit_size);
    }

    template <typename T>
    proc pow_trivial(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t power, Optims const& optims) -> Slice<T>
    {
        assert(is_striped_number(num));

        assert(are_aliasing<T>(*to, num) == false);
        assert(are_aliasing<T>(*aux, num) == false);
        assert(are_aliasing<T>(*aux, *to) == false);

        const size_t required_to_size = required_pow_out_size(log2(num), power, BIT_SIZE<T>);
        const size_t required_aux_size = required_pow_trivial_aux_size(log2(num), power, BIT_SIZE<T>);
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
    proc pow(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t power, Optims const& optims) -> Slice<T>
    {
        size_t bitsize = log2(num);
        size_t at_least_to = required_pow_out_size(bitsize, power, BIT_SIZE<T>);
        size_t at_least_aux = required_pow_trivial_aux_size(bitsize, power, BIT_SIZE<T>);
        assert(to->size >= at_least_to);
        assert(aux->size >= at_least_aux);

        if(power < optims.pow_trivial_below_power)
            return pow_trivial<T>(to, aux, num, power, optims);

        const size_t required_aux_size = required_pow_by_squaring_aux_size(bitsize, power, BIT_SIZE<T>);

        if(aux->size >= required_aux_size)
            return pow_by_squaring(to, aux, num, power, optims);

        return pow_trivial<T>(to, aux, num, power, optims);
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

    template <typename T>
    func single_root_shifting(T x, T n) -> T
    {
        if(n == 0)
            return 1;
        if (x <= 1) 
            return x;

        T r = 1;
        T s = ((find_last_set_bit(x) / n) * n);

        while(true)
        {
            if(s < n)
                break;

            s -= n; 
            r <<= 1; 
            T power = single_pow_by_squaring<T>(r | 1, n);
            T bit = power <= (x >> s);
            r |= bit;
        }

        return r;
    }

    template <typename T>
    func single_root_newton(T of_value, T power) -> T
    {
        const T x = of_value; 
        const T n = power;

        if(x == 0)
            return 1;
        if (x <= 1) 
            return x;

        //T u = x;
        //T s = x+1;

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

        const T log2x = find_last_set_bit(x);

        const T upper_initial_estimate = cast(T) 1 << ((log2x + n) / n);
        T r = upper_initial_estimate;
        T prev_r = -1;

        while(true)
        {
            // the newton update
            T new_upper_r = (n-1) * r + x / single_pow_by_squaring(r, n-1);

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
    proc root(Slice<T>* to, Slice<T>* aux, Slice<const T> num, size_t root, Optims const& optims) -> Slice<T>
    {
        size_t num_bit_size = log2(num);   
        size_t required_to_size = required_root_out_size(num.size, root, BIT_SIZE<T>);
        size_t required_aux_size = required_root_aux_size(num_bit_size, root, BIT_SIZE<T>);

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
            to->data[0] = cast(T) single_root_shifting(num[0], cast(T) root);
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

            //new_upper_r = (root-1) * prev_estimate + num / single_pow_by_squaring(prev_estimate, root-1);

            
            //the sizes we will need to store the results:
            size_t powed_size = required_pow_out_size(log2(prev_estimate), root - 1, BIT_SIZE<T>);
            size_t muled_size = required_mul_out_size(prev_estimate.size, 1);

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


            Div_Result<T> res = div<T>(&div_to, &rem_to, num, powed, optims);
            Slice<T> dived = res.quotient;
            assert(is_striped_number(dived));

            size_t add_size = required_add_out_size(muled_size, dived.size);

            Slice<T> new_estimate_to = out1;
            assert(new_estimate_to.data != prev_estimate.data && "storages must differ");
            assert(new_estimate_to.size >= add_size || true &&  
                "should fit into the prepared slot - maybe will fail because add_size is theoretical reuquired size"
                "but new_estimate_to.size is the actual real size which might - and probably be - smaller than the theoretical size"
                "we might have to lie about size to the algorhitm here!");

            //we perform the multiply add fused
            Batch_Overflow<T> fused_result;
            
            if(dived.size != 0)
                fused_result = batch_fused_mul_add_overflow<T>(&new_estimate_to, dived, prev_estimate, cast(T) root - 1, optims);
            else
                fused_result = batch_mul_overflow<T>(&new_estimate_to, prev_estimate, cast(T) root - 1, optims);
                
            assert(fused_result.overflow == 0 && "should not overflow");
            assert(is_striped_number(fused_result.slice));
            std::swap(out1, out2);


            //r = new_upper_r / n;
            Slice<T> new_upper_estimate = fused_result.slice;
            Slice<T> new_estimate_pre = batch_div_overflow<T>(&new_upper_estimate, new_upper_estimate, cast(T) root, optims).slice;
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
    func single_log_heyley(T of_value, T base) -> T
    {
        const T x = of_value;
        const T b = base;

        //b^y = x
        //y = log_b(x)

        if(x <= 1)
            return 0;
        if(base <= 1)
            return 0;

        const T log2x = find_last_set_bit(x);
        const T log2b = find_last_set_bit(b);

        assert(log2x != 0);
        assert(log2b != 0);

        const T lower_initial_estimate = log2x / (log2b + 1);
        const T higher_initial_estimate = (log2x + log2b) / (log2b);
        T y = higher_initial_estimate;
        T prev_y = y;

        while(true)
        {
            T value_from_curr_approximate = single_pow_by_squaring(b, y);
            T c_x = value_from_curr_approximate;

            prev_y = y;
            if(c_x <= x)
                break;

            T new_y = y - 2*(c_x - x) / (x + c_x);

            //if we didnt move at all we are one above the answer
            if(new_y == y)
                return new_y - 1;

            y = new_y;
        }

        return prev_y;

    }

    template <typename T>
    func single_log_bsearch(T of_value, T base) -> T
    {
        const T x = of_value;
        const T b = base;

        if(x <= 1)
            return 0;
        if(base <= 1)
            return 0;

        const T log2x = find_last_set_bit(x);
        const T log2b = find_last_set_bit(b);

        assert(log2x != 0);
        assert(log2b != 0);

        const T lower_initial_estimate = log2x / (log2b + 1);
        const T higher_initial_estimate = (log2x + log2b) / (log2b);

        T curr_lo = lower_initial_estimate;
        T curr_hi = higher_initial_estimate;

        while(true)
        {
            T mid = (curr_lo + curr_hi) / 2;
            T curr_approx = single_pow_by_squaring(b, mid);

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
    proc reverse(Slice<T>* arr) -> void
    {
        Slice<T> ref = *arr;
        size_t half_size = ref.size / 2;
        for(size_t i = 0; i < half_size; i++)
            std::swap(ref[i], ref[ref.size - i - 1]);
    }

    struct Power_And_Powed
    {
        size_t power = 0;
        size_t powed = 0;

        bool constexpr operator ==(Power_And_Powed const&) const noexcept = default;
    };

    func highest_power_of_in_type(size_t type_bit_size, size_t base) -> Power_And_Powed
    {
        assert(type_bit_size <= 64 && "no higher sizes supported");

        if(type_bit_size <= 1) 
            return {0, 1};

        switch(base)
        {
        case 10: 
            if(type_bit_size < 4)    return {0, 1};
            if(type_bit_size == 4)   return {1, 10};
            if(type_bit_size == 8)   return {2, 100};
            if(type_bit_size == 16)  return {4, 10'000};
            if(type_bit_size == 32)  return {9, 1'000'000'000};
            if(type_bit_size == 64)  return {19, 10'000'000'000'000'000'000};

            break;

        case 2:
            if(type_bit_size == 64)
                type_bit_size = 63;

            return {(type_bit_size - 1), 1ull << (type_bit_size - 1)};

        case 16:
            return {(type_bit_size - 1)/4, 1ull << (((type_bit_size - 1)/4) * 4)};
        }

        size_t max = 0;
        if(type_bit_size < 64)
            max = (1ull << type_bit_size) - 1;
        else
            max = cast(size_t) -1;

        size_t power = single_log_bsearch<size_t>(max, base);
        size_t powed = single_pow_by_squaring<size_t>(base, power);

        return {power, powed};
    }

    func required_to_base_out_size(size_t num_size, size_t to_base, size_t single_digit_bit_size) -> size_t
    {
        Power_And_Powed highest = highest_power_of_in_type(single_digit_bit_size / 2, to_base); 
        
        //into a single digit fits `power` powers of to_base
        // => we need to `power + 1` to_base to represent any number of the digit
        // => into the whole number fit `num_size * power` powers
        // => we need `num_size * (power + 1)` digits
        //@TODO: why do we need the +highest.power?
        return num_size * 2 * (highest.power + 1) + highest.power;
    }

    auto required_from_base_out_size(size_t rep_size, size_t from_base, size_t single_digit_bit_size) -> size_t
    {
        Power_And_Powed highest = highest_power_of_in_type(single_digit_bit_size / 2, from_base); 

        //[][][][][]
        //----|----| 
        //  => 3 rep digits per num digit
        //  => we just need num digits / rep digits per digit rounded up
        return div_round_up(rep_size, highest.power);
    }

    template <typename Num>
    struct Optional_Slice
    {
        bool ok;
        Slice<Num> slice;
    };

    template <typename Num>
    struct Pop_Digit
    {
        Slice<Num> buffer;
        Slice<const Num> number; //can also point into a different buffer from buffer
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

    template <typename Num>
    func pop_digit_init(Slice<Num>* buffer, Slice<const Num> num, Num base) -> Pop_Digit<Num>
    {
        Pop_Digit popper = Pop_Digit{*buffer, num, cast(size_t) base};

        popper.highest = highest_power_of_in_type(HALF_BIT_SIZE<Num>, cast(size_t) base);
        if(high_bits(base) != 0)
        {
            popper.error = Pop_Digit<Num>::BASE_TOO_LARGE;
            popper.ok = false;
        }

        if(base < 2)
        {
            popper.error = Pop_Digit<Num>::BASE_TOO_SMALL;
            popper.ok = false;
        }

        if(buffer->size < num.size)
        {
            popper.error = Pop_Digit<Num>::SMALL_BUFFER;
            popper.ok = false;
        }

        if(is_striped_number(num) == false)
        {
            popper.error = Pop_Digit<Num>::NON_STRIPPED_NUM;
            popper.ok = false;
        }

        return popper;
    }

    template <typename Num>
    proc pop_digit(Pop_Digit<Num>* state, Optims const& optims) -> Num
    {
        static_assert(is_allowed_digit<Num>, "!");

        if(state->stashed_power == 0)
        {
            if(state->number.size == 0)
            {
                state->ok = false; 
                return -1;
            }

            Batch_Overflow<Num> res = batch_div_overflow<Num>(&state->buffer, state->number, cast(Num) state->highest.powed, optims);

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

        Num digit = cast(Num) state->stashed_value % state->stashed_base;
        state->stashed_value /= state->stashed_base;
        state->stashed_power -= 1;

        return digit;
    }

    template <typename Num>
    struct Push_Digit
    {
        Slice<Num> buffer;
        size_t num_size = 0;
        size_t stashed_base = 0;
        size_t stashed_digit = 0;
        size_t stashed_power = 0;
        Power_And_Powed highest;
    };

    template <typename Num>
    func push_digit_init(Slice<Num>* buffer, size_t num_size = 0) -> Push_Digit<Num>
    {
        Push_Digit pusher = Push_Digit{*buffer, num_size};
        return pusher;
    }

    template <typename Num>
    func push_digit(Push_Digit<Num>* state, Num digit, Num base, Optims const& optims, bool flush = false) -> bool
    {
        static_assert(is_allowed_digit<Num>, "!");
        if(base < 2 || digit > base || state->num_size + 1 > state->buffer.size)
            return false;
        
        if(base != state->stashed_base)
        {
            if(state->stashed_power != 0)
                if(push_digit<Num>(state, 0, cast(Num) state->stashed_base, optims, true))
                    return false;

            state->highest = highest_power_of_in_type(BIT_SIZE<Num>, base);
            state->stashed_base = cast(size_t) base;
        }

        if(state->num_size + 1 > state->buffer.size)
            assert(false && "shouldnt be possible");

        //flush only attempts to flush stashed_digit into the num buffer
        if(flush == false)
        {
            state->stashed_digit *= base;
            state->stashed_digit += cast(size_t) digit;
            state->stashed_power += 1;

            if(state->stashed_power < state->highest.power)
                return state;
        }

        Num powed = cast(Num) state->highest.powed;
        if(flush)
            powed = single_pow_trivial<Num>(base, cast(Num) state->stashed_power);

        Slice<Num> curr_num = trim(state->buffer, state->num_size);
        size_t overflown_times = 0;

        let mul_res = batch_mul_overflow<Num>(&curr_num, curr_num, powed, optims);
        if(mul_res.overflow != 0)
        {
            overflown_times += 1;
            state->num_size += 1;
            curr_num = trim(state->buffer, state->num_size);
            curr_num[state->num_size - 1] = mul_res.overflow;
        }

        let add_res = batch_add_overflow_short<Num>(&curr_num, curr_num, cast(Num) state->stashed_digit, Location::IN_PLACE);
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

    template <typename Num>
    proc push_digit_result(Push_Digit<Num>* state, Optims const& optims) -> Optional_Slice<Num>
    {
        if(state->stashed_power != 0)
        {
            if(push_digit<Num>(state, 0, cast(Num) state->stashed_base, optims, true) == false)
                return Optional_Slice<Num>{false};
        }

        return Optional_Slice<Num>{true, trim(state->buffer, state->num_size)};
    }
    
    template <typename Num, typename Rep, typename Conversion_Fn>
    proc from_base(Slice<Num>* num, Slice<const Rep> rep, Num base, Conversion_Fn conversion, Optims const& optims) -> Optional_Slice<Num>
    {
        static_assert(is_allowed_digit<Num>, "!");
        static_assert(sizeof(Num) >= sizeof(Rep), "!");
        if constexpr(std::is_same_v<Rep, Num>)
            assert(check_are_aliasing<Rep>(rep, *num) == false);

        assert(required_from_base_out_size(rep.size, base, BIT_SIZE<Num>) <= num->size 
            && "there must be enough size in num to represent the number even in the worst case scenario");

        mut state = push_digit_init(num);
        
        for(size_t i = 0; i < rep.size; i++)
        {
            Num digit = conversion(rep[i]);
            if(push_digit(&state, digit, base, optims) == false)
                return {false};
        }

        return push_digit_result(&state, optims);
    }

    template <typename Num, typename Rep, typename Conversion_Fn>
    proc to_base_ensembled(Slice<Rep>* rep, Slice<Num>* temp, Slice<const Num> num, Num base, Conversion_Fn conversion, Optims const& optims) -> Slice<Rep>
    {
        static_assert(is_allowed_digit<Num>, "!");
        if constexpr(std::is_same_v<Num, Rep>)
            assert(check_are_aliasing<Num>(num, *rep) == false);

        if constexpr (DO_RUNTIME_ONLY)
        {
            size_t at_least = required_to_base_out_size<Num>(num.size, base);
            assert(at_least <= rep->size && "there must be enough size rep represent the number even in the worst case scenario");
        }

        mut state = pop_digit_init<Num>(temp, num, base);
        size_t i = 0;
        for(; ; i++)
        {
            Num popped = pop_digit(&state, optims);
            if(state.ok == false)
                break;

            Rep converted = conversion(popped);
            (*rep)[i] = converted;
        }

        assert(state.error == Pop_Digit<Num>::OK);
        Slice<Rep> converted = trim(*rep, i);
        reverse(&converted);
        return converted;
    }

    template <typename Num, typename Rep, typename Conversion_Fn>
    proc to_base(Slice<Rep>* rep, Slice<Num>* temp, Slice<const Num> num, Num base, Conversion_Fn conversion, Optims const& optims) -> Slice<Rep>
    {
        static_assert(is_allowed_digit<Num>, "!");
        static_assert(sizeof(Rep) <= sizeof(Num), "!");

        if constexpr(std::is_same_v<Num, Rep>)
            assert(check_are_aliasing<Num>(num, *rep) == false);

        assert(is_striped_number(num));
        assert(temp->size >= num.size && "temp must be big enough");
        assert(base >= 2 && "base must be bigger than two");
        assert(high_bits(base) == 0 && "base must be low bits only of Num type (so that short div algorhimt can be used)");

        size_t at_least = required_to_base_out_size(num.size, base, BIT_SIZE<Num>);
        assert(at_least <= rep->size && "there must be enough size rep represent the number even in the worst case scenario");

        bool is_first = true;
        size_t to_size = 0;

        if(num.size == 0)
            return trim(*rep, 0);


        Power_And_Powed highest = highest_power_of_in_type(HALF_BIT_SIZE<Num>, base);
        assert(highest.powed > 1);

        Slice<const Num> val_current = num;
        bool iterate = true;
        while(iterate)
        {
            Batch_Overflow<Num> res = batch_div_overflow<Num>(temp, val_current, cast(Num) highest.powed, optims);

            Num highest_power_rem = res.overflow;
            val_current = striped_trailing_zeros(res.slice);

            bool is_last_block = val_current.size == 0;
            for(size_t j = 0; j < highest.power; j++)
            {
                if(is_last_block && highest_power_rem == 0)
                {
                    iterate = false;
                    break;
                }

                Num rem = highest_power_rem % base;
                Rep digit = conversion(rem);
                (*rep)[to_size] = digit;
                to_size++;

                highest_power_rem /= base;
            }


            is_first = false;
        }

        Slice<Rep> converted = trim(*rep, to_size);
        reverse(&converted);
        return converted;
    }

    func single_factorial(size_t of_value) -> size_t
    {
        size_t calculated = 1;
        for(size_t factor = 2; factor <= of_value; factor++)
            calculated *= factor;

        return calculated;
    }

    static constexpr size_t MAX_NATIVE_FACTORIAL = 20;

    func single_factorial_fast(size_t of_value) -> size_t
    {
        switch(of_value) {
            case 0: return single_factorial(0);
            case 1: return single_factorial(1);
            case 2: return single_factorial(2);
            case 3: return single_factorial(3);
            case 4: return single_factorial(4);
            case 5: return single_factorial(5);
            case 6: return single_factorial(6);
            case 7: return single_factorial(7);
            case 8: return single_factorial(8);
            case 9: return single_factorial(9);
            case 10: return single_factorial(10);
            case 11: return single_factorial(11);
            case 12: return single_factorial(12);
            case 13: return single_factorial(13);
            case 14: return single_factorial(14);
            case 15: return single_factorial(15);
            case 16: return single_factorial(16);
            case 17: return single_factorial(17);
            case 18: return single_factorial(18);
            case 19: return single_factorial(19);
            case 20: return single_factorial(20);
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

    template <typename T>
    proc factorial(Slice<T>* out, size_t of_value, Optims const& optims) -> Slice<T>
    {
        assert(out->size >= required_factorial_size(of_value));

        bool did_finish = true;
        if(of_value < MAX_NATIVE_FACTORIAL)
        {
            size_t fact = single_factorial_fast(of_value);
            return from_number(out, fact);
        }

        constexpr size_t highest_possible_native = single_factorial(MAX_NATIVE_FACTORIAL);
        Slice<T> current_value = from_number(out, highest_possible_native);
        for(size_t factor = MAX_NATIVE_FACTORIAL + 1; factor <= of_value; factor++)
            current_value = mul<T>(out, current_value, cast(T) factor, optims);

        return current_value;
    }

    #undef let
    #undef mut
    #undef proc
    #undef func
    #undef cast
}