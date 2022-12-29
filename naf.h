#pragma once

#include <algorithm>
#include <cassert>
#include <exception>
#include <execution>
#include <iostream>
#include <ranges>
#include <vector>
#include <bit>

//define c restruct keyword for different compilers
// if not supported compiler just use nothing
#if defined(__clang__)
    #define restrict __restrict__
#elif defined(__GNUC__)
    #define restrict __restrict__
#elif defined(_MSC_VER)
    #define restrict __restrict
#else
    #define restrict
#endif

namespace naf
{
    constexpr size_t SIMD_SIZE = 256; //SSE 

    constexpr size_t div_round_up(size_t nom, size_t den) noexcept
    {
        return (nom + den - 1) / den;
    }

    template <typename T, typename Alloc = std::allocator<T>>
    struct Vector : Alloc
    {
        static_assert(std::is_scalar_v<T> || (std::is_trivially_copyable_v<T> && std::is_trivially_constructible_v<T>), "must be trivial type!");

        T* _data = nullptr;
        size_t _size = 0;
        size_t _capacity = 0;

        static constexpr size_t BLOCK_SIZE = SIMD_SIZE / sizeof(T); 

        using value_type      = T;
        using size_type       = size_t;
        using difference_type = ptrdiff_t;
        using pointer         = T*;
        using const_pointer   = const T*;
        using reference       = T&;
        using const_reference = const T&;
        using iterator_category = std::contiguous_iterator_tag;

        using iterator       = T*;
        using const_iterator = const T*;

        using reverse_iterator       = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

        constexpr auto begin() noexcept -> iterator {
            return this->_data;
        }

        constexpr auto begin() const noexcept -> const_iterator {
            return this->_data;
        }

        constexpr auto end() noexcept -> iterator{
            return this->_data + this->_size;
        }

        constexpr auto end() const noexcept -> const_iterator {
            return this->_data + this->_size;
        }

        constexpr auto operator[](size_t index) const noexcept -> T const& { 
            assert(0 <= index && index < this->_size && "index out of range"); 
            return this->_data[index]; 
        }
        constexpr auto operator[](size_t index) noexcept -> T& { 
            assert(0 <= index && index < this->_size && "index out of range"); 
            return this->_data[index]; 
        }

        constexpr auto data() noexcept -> iterator{
            return this->_data;
        }

        constexpr auto data() const noexcept -> const_iterator {
            return this->_data;
        }

        constexpr auto size() const noexcept -> size_t {
            return this->_size;
        }

        constexpr auto capacity() const noexcept -> size_t{
            return this->_capacity;
        }


        constexpr auto back() const noexcept -> T {
            assert(_size > 0);
            return _data[_size - 1];
        }

        constexpr auto front() const noexcept -> T{
            assert(_size > 0);
            return _data[0];
        }

        constexpr void swap(Vector& other) noexcept
        {
            std::swap(this->_data, other._data);
            std::swap(this->_size, other._size);
            std::swap(this->_capacity, other._capacity);
        }
        constexpr auto allocator() -> Alloc* {
            return (Alloc*) this;
        
        }

        constexpr Vector() noexcept = default;
        constexpr Vector(Vector&& other) noexcept
        {
            swap(other);
        }


        constexpr Vector(Vector const& other)
        {
            if(other._size == 0)
                return;

            _copy_other_with_capacity(other, other._size);
        }

        constexpr Vector& operator=(Vector&& other) noexcept
        {
            swap(other);
            return *this;
        }

        constexpr Vector& operator=(Vector const& other)
        {
            if(&other == this)
                return *this;

            resize_for_overwrite(other._size);
            memcpy(this->_data, other._data, other._size*sizeof(T));
            return *this;
        }

        constexpr ~Vector() noexcept
        {
            if(this->_capacity == 0)
                return;

            allocator()->deallocate(this->_data, this->_capacity);
        }

        constexpr void _copy_other_with_capacity(Vector const& other, size_t capacity)
        {
            assert(capacity >= other._size);
            this->_capacity = capacity;
            this->_size = other._size;
            this->_data = allocator()->allocate(this->_capacity);

            memcpy(this->_data, other._data, other._size*sizeof(T));
        }

        constexpr void reserve(size_t to_capacity)
        {
            if(to_capacity <= this->_capacity)
                return;

            size_t ceiled = std::bit_ceil(to_capacity);
            size_t rounded = div_round_up(ceiled, BLOCK_SIZE)*BLOCK_SIZE;

            Vector copy;
            copy._copy_other_with_capacity(*this, rounded);
            swap(copy);
        }

        constexpr void resize_for_overwrite(size_t to_size)
        {
            reserve(to_size);
            this->_size = to_size;
        }

        constexpr void resize(size_t to_size)
        {
            reserve(to_size);
            
            if(to_size > this->_size)
                memset(this->_data + this->_size, 0, (to_size - this->_size)*sizeof(T));

            this->_size = to_size;
        }

        constexpr void push_back(T what)
        {
            reserve(this->_capacity + 1);
            this->_data[this->_size] = what;
            this->_size ++;
        }

        constexpr void clear() noexcept
        {
            this->_size = 0;
        }

        constexpr void pop_back() noexcept
        {
            assert(_size > 0);
            this->_size --;
        }
    };
}
namespace naf
{
    using std::ranges::iota_view;
    using std::views::reverse;

    using Digit = std::int64_t;
    using Digits = Vector<Digit>;
    using Size = std::size_t;

    const auto parallel = std::execution::par;

    /// The radix used for the integers.
    constexpr static const Digit base = Digit(1) << 31;
    constexpr static const Digit half_base = Digit(1) << 30;
    /// The individual digits can range from -max_digit to max_digit inclusive.
    constexpr static const Digit max_digit = base - 1;

    namespace ops
    {
        static Digit add(Digit x, Digit y) { return x + y; };

        static Digit sub(Digit x, Digit y) { return x - y; }; 
        static Digit sub_reverse(Digit x, Digit y) { return y - x; }; 

        static Digit add_zero(Digit y) { return 0 + y; };
        static Digit sub_zero(Digit y) { return 0 - y; };
    }

    static Size find_first_non_zero(Digits const& digits) noexcept
    {
        //@TODO make paralel here
        Size i = digits.size();
        for(; i-- > 0;)
            if(digits[i] != 0)
                break;

        return i;
    }

    //we could probably come up with a procedure that merges both inc_or_decr and add_or_sub but 
    // the code woiuld just be more complex and would reuqire us passing iterators

    /// Uses Avizienis' parallel addition algorithm, as described on page 4 in
    /// https://arxiv.org/pdf/1102.5683.pdf
    static void incr_or_decr(Digits* restrict carry, Digits* restrict left, Digits const& restrict right, bool is_increment)
    {
        const auto op_binary = is_increment ? ops::add : ops::sub;
        Size max_size = std::max(left->size(), right.size());

        //left->resize(max_size + 1);
        left->resize(max_size);


        std::transform(parallel, left->begin(), left->begin() + right.size(),
            right.begin(), left->begin(), op_binary);

        //@BUG: ?????
        //carry->clear();
        //carry->resize(max_size);
        //std::transform(parallel, 
            //left->begin() + 1, left->begin() + 1 + carry->size(), 
            //carry->begin(), left->begin() + 1, op_binary);

        Size non_zero = find_first_non_zero(*left);
        if(non_zero + 1 < left->size())
            left->resize(non_zero + 1);
    }

    static void add_or_sub(Digits* restrict output, Digits* restrict carry, Digits const& restrict left, Digits const& restrict right, bool is_addition)
    {
        const auto op_binary = is_addition ? ops::add : ops::sub;
        const auto op_unary = is_addition ? ops::add_zero : ops::sub_zero;

        Size max_size = std::max(left.size(), right.size());
        Size min_size = std::min(left.size(), right.size());

        //output->clear();
        //output->resize(max_size + 1);
        output->resize_for_overwrite(max_size);

        //perform the binary op on common portion of all inputs
        std::transform(parallel, left.begin(), left.begin() + min_size,
            right.begin(), output->begin(), op_binary);

        //patch up the remainign size:
        //if left operand remains just copy the elements over
        if(left.size() > min_size)
        {
            std::copy(parallel, left.begin() + min_size, left.end(), 
                output->begin() + min_size);
        }
        //if right operand remains perform transform the elements using unary op
        else
        {
            std::transform(parallel, right.begin() + min_size, right.end(), 
                output->begin() + min_size, op_unary);
        }

        //@BUG: ?????
        //carry->clear();
        //carry->resize(max_size);
        //std::transform(parallel, 
            //output->begin() + 1, output->begin() + 1 + carry->size(), 
            //carry->begin(), output->begin() + 1, op_binary);

        Size non_zero = find_first_non_zero(*output);
        if(non_zero + 1 < output->size())
            output->resize(non_zero + 1);
    }

    /// Calculate the additive inverse of an Integer in place.
    static void negate(Digits* digits) noexcept {
        std::transform(parallel, digits->begin(), digits->end(), digits->begin(),
            [](Digit digit) { return -digit; });
    }

    static void negate(Digits* restrict output, Digits const& restrict from) {
        //output->clear();
        //output->resize(from.size());
        output->resize_for_overwrite(from.size());
        std::transform(parallel, from.begin(), from.end(), output->begin(),
            [](Digit digit) { return -digit; });
    }

    static Digit compare_zero(Digits const& digits) noexcept
    {
        if(digits.size() == 0)
            return 0;

        return digits.back();
    }

    static Digit compare(Digits* restrict temp, Digits* restrict carry, Digits const& restrict left, Digits const& restrict right) 
    {
        add_or_sub(temp, carry, left, right, false);
        return compare_zero(*temp);
    }

    static bool are_equal(Digits const& restrict left, Digits const& restrict right) noexcept
    {
        if(left.size() != right.size())
            return false;

        return memcmp(left.data(), right.data(), right.size() * sizeof(Digit)) == 0;
    }

    static Digit compare(Digits const& restrict left, Digits const& restrict right) 
    {
        Digits carry; 
        Digits temp; 
        return compare(&temp, &carry, left, right);
    }

    static void mul( 
        Digits* restrict output, Digits* restrict temp,
        Digits* restrict augend, Digits* restrict carry, 
        Digits const& restrict left, Digits const& restrict right)
    {
        // Optimization: this algorithm is faster
        // when multiplying a shorter number by a longer number.
        if (left.size() > right.size()) 
            return mul(output, temp, augend, carry, right, left);

        Size result_length = left.size() + right.size();

        //augend->clear();
        //carry->clear();
        augend->resize_for_overwrite(result_length);
        carry->resize_for_overwrite(result_length);
        output->clear();
        output->reserve(result_length);

        const auto this_indices = iota_view(Size(0), left.size());
        const auto other_indices = iota_view(Size(0), right.size());

        std::for_each(this_indices.begin(), this_indices.end(), [&](Size i) {
            std::fill(parallel, augend->begin(), augend->end(), 0);
            std::fill(parallel, carry->begin(), carry->end(), 0);

            //can be paralel but requires random access iota_view @TODO
            std::for_each(other_indices.begin(), other_indices.end(),
                [&](Size j) {
                    Digit product = left[i] * right[j];
                    (*augend)[i + j] = product % naf::base;
                    (*carry)[i + j + 1] = product / naf::base;
                });

            naf::incr_or_decr(temp, output, *augend, true);
            naf::incr_or_decr(temp, output, *carry, true);
        });
    }

    static void half(Digits* restrict output, Digits const& restrict input) 
    {
        //output->clear();
        //output->resize(input.size());
        output->resize_for_overwrite(input.size());

        if(input.size() == 0)
            return;

        std::transform(parallel, input.begin(), input.end(), output->begin(),
            [](Digit digit) { return digit / 2; });

        std::transform(parallel, input.begin() + 1, input.end(), output->begin(),
            output->begin(), [](Digit high, Digit low) {
                return low + (high % 2) * half_base;
            });

        if (output->back() == 0) {
            output->pop_back();
        }
    }
}