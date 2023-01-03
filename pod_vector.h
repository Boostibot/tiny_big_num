#pragma once

#include <memory>

template <typename T, typename Alloc = std::allocator<T>>
struct POD_Vector : Alloc
{
    static_assert(std::is_scalar_v<T> || (std::is_trivially_copyable_v<T> && std::is_trivially_constructible_v<T>), "must be trivial type!");

    T* _data = nullptr;
    size_t _size = 0;
    size_t _capacity = 0;

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

    constexpr auto empty() const noexcept -> bool{
        return _size == 0;
    }

    constexpr void swap(POD_Vector& other) noexcept
    {
        std::swap(this->_data, other._data);
        std::swap(this->_size, other._size);
        std::swap(this->_capacity, other._capacity);
    }
    constexpr auto allocator() -> Alloc* {
        return (Alloc*) this;

    }

    constexpr POD_Vector() noexcept = default;
    constexpr POD_Vector(POD_Vector&& other) noexcept
    {
        swap(other);
    }

    constexpr POD_Vector(POD_Vector const& other)
    {
        if(other._size == 0)
            return;

        _copy_other_with_capacity(other, other._size);
    }

    constexpr POD_Vector(size_t size)
    {
        resize(size);
    }

    constexpr POD_Vector& operator=(POD_Vector&& other) noexcept
    {
        swap(other);
        return *this;
    }

    constexpr POD_Vector& operator=(POD_Vector const& other)
    {
        if(&other == this)
            return *this;

        resize_for_overwrite(other._size);
        memcpy(this->_data, other._data, other._size*sizeof(T));
        return *this;
    }

    constexpr ~POD_Vector() noexcept
    {
        if(this->_capacity == 0 || this->_data == nullptr)
            return;

        allocator()->deallocate(this->_data, this->_capacity);
    }

    constexpr void _copy_other_with_capacity(POD_Vector const& other, size_t capacity)
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
        ceiled = (ceiled < 8) ? 8 : ceiled;

        POD_Vector copy;
        copy._copy_other_with_capacity(*this, ceiled);
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
        reserve(this->_size + 1);
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
