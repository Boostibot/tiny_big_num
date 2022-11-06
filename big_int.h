#pragma once

#include "preface.h"

template <typename Grow>
concept grower = requires(size_t to_fit, size_t capacity, size_t size, size_t static_capacity, size_t elem_size)
{
    { Grow::run(to_fit, capacity, size, static_capacity, elem_size) } -> std::convertible_to<size_t>;
};

template<size_t mult_, size_t add_, size_t base_elems_>
struct Def_Grow
{
    static constexpr size_t mult = mult_;
    static constexpr size_t add = add_;
    static constexpr size_t base_elems = base_elems_;

    static constexpr size_t run(size_t to_fit, size_t capacity, size_t size, size_t static_capacity, size_t elem_size)
    {
        cast(void) size;
        cast(void) static_capacity;
        cast(void) elem_size;

        size_t realloc_to = capacity 
            ? capacity * mult + add 
            : base_elems;

        while(realloc_to < to_fit)
            realloc_to = realloc_to * mult + add;

        return realloc_to;
    }
};

static_assert(grower<Def_Grow<2, 0, 8>>);

namespace detail
{   
    template<integral T, size_t static_capacity_>
    struct Big_Int_Data
    {
        T* data = nullptr;
        size_t size = 0;
        size_t capacity = 0;
        char sign = false; //positive
                           //will maybe get fused with capacity as last bit (to revove padding and enable more elemensts to 
                           // be stored inline inside the struct itself
        T static_data_[static_capacity_] = {0};
    };

    template<integral T>
    struct Big_Int_Data<T, 0>
    {
        T* data = nullptr;
        size_t size = 0;
        size_t capacity = 0;
        char sign = false; //positive
    };

}

template <
    integral T, 
    size_t static_capacity_, 
    allocator Alloc_, 
    grower Grow
>
struct Big_Int_;

#define BIG_INT_TEMPL class T, size_t scap, class Alloc, class Grow
#define Big_Int_T Big_Int_<T, scap, Alloc, Grow>

template<BIG_INT_TEMPL>
func static_data(Big_Int_T* data) -> T*;

template<BIG_INT_TEMPL>
func static_data(const Big_Int_T* data) -> const T*;

namespace detail 
{
    template<BIG_INT_TEMPL>
    proc assign(Big_Int_T* to, Big_Int_T const& other) -> void;

    template<BIG_INT_TEMPL>
    proc swap(Big_Int_T* left, Big_Int_T* right) noexcept -> void;

    template<BIG_INT_TEMPL>
    proc alloc_data(Big_Int_T* num, size_t capacity) -> void;

    template<BIG_INT_TEMPL>
    proc dealloc_data(Big_Int_T* num) -> void;

    template<BIG_INT_TEMPL, typename T_>
    proc to_owned(Big_Int_T* num, span<T_> const& other) -> void;
}



using Max_Unsigned_Type = u64;

template <
    integral T_, 
    size_t static_capacity_ = 2,
    allocator Alloc_ = std::allocator<T_>, 
    grower Grow = Def_Grow<2, 0, 8>
>
struct Big_Int_ : Alloc_, detail::Big_Int_Data<T_, static_capacity_>
{   
    using T = T_;
    using Alloc = Alloc_;
    using Big_Int_Data     = detail::Big_Int_Data<T, static_capacity_>;

    using allocator_type   = Alloc;
    using slice_type       = span<T>;
    using const_slice_type = span<const T>;
    using grow_type        = Grow;

    constexpr static size_t static_capacity = cast(size_t) static_capacity_; 
    constexpr static bool has_static_storage = static_capacity != 0;
    constexpr static bool has_stateful_alloc = !std::is_empty_v<Alloc>;

    constexpr Big_Int_() = default;
    constexpr Big_Int_(Big_Int_&& other, Alloc alloc = Alloc()) noexcept
        : Alloc{move(alloc)} { detail::swap(this, &other); }

    constexpr Big_Int_(Max_Unsigned_Type val, Alloc alloc = Alloc()) 
        : Alloc{move(alloc)} 
    { 
        constexpr size_t count = sizeof(Max_Unsigned_Type) / sizeof(T);

        detail::alloc_data(this, count);

        size_t i = 0;
        for(; i < count; i++)
        {
            Max_Unsigned_Type curr_digit = val >> (i * sizeof(T) * CHAR_BIT);
            if(curr_digit == 0)
                break;

            this->data[i] = cast(T) curr_digit;
        }

        this->size = i;
    }

    explicit constexpr Big_Int_(slice_type slice, Alloc alloc = Alloc()) 
        : Alloc{move(alloc)} { detail::to_owned(this, slice); }

    explicit constexpr Big_Int_(const_slice_type slice, Alloc alloc = Alloc()) 
        : Alloc{move(alloc)} { detail::to_owned(this, slice); }

    constexpr Big_Int_(T* data, size_t size, size_t capacity, Alloc alloc, const T (&static_data)[static_capacity]) noexcept 
        requires (has_static_storage)
    : Alloc{move(alloc)}, Big_Int_Data{data, size, capacity} 
    {
        copy_n(static_data(this), static_data, static_capacity);
    }

    constexpr Big_Int_(T* data, size_t size, size_t capacity, Alloc alloc = Alloc()) noexcept
        : Alloc{move(alloc)}, Big_Int_Data{data, size, capacity} {}

    constexpr Big_Int_(Alloc alloc) noexcept
        : Alloc{move(alloc)} {}

    constexpr Big_Int_(const Big_Int_& other)
        requires std::is_copy_constructible_v<Alloc>
    : Alloc{cast(Alloc) other}
    {
        detail::to_owned(this, cast(slice_type) other);
    }

    constexpr ~Big_Int_() noexcept
    {
        detail::dealloc_data(this);
    }

    constexpr Big_Int_& operator=(const Big_Int_& other)
    {
        detail::assign(this, other);
        return *this;
    }

    constexpr Big_Int_& operator=(Big_Int_&& other) noexcept
    {
        detail::swap(this, &other);
        return *this;
    }

    constexpr operator slice_type() noexcept 
    { 
        return slice_type{this->data, this->size}; 
    }

    constexpr operator const_slice_type() const noexcept 
    { 
        return const_slice_type{this->data, this->size}; 
    }

    func operator[](size_t i) noexcept -> T& 
    {
        assert(i < this->size);
        return this->data[i];
    }

    func operator[](size_t i) const noexcept -> T const& 
    {
        assert(i < this->size);
        return this->data[i];
    }
};

namespace std 
{
    template <BIG_INT_TEMPL>
    proc swap(Big_Int_T& num1, Big_Int_T& num2)
    {
        return detail::swap(&num1, &num2);
    }
}

template <BIG_INT_TEMPL>
func static_data(Big_Int_T* data) -> T*
{
    if constexpr(scap != 0)
        return data->static_data_;
    else
        return nullptr;
}

template <BIG_INT_TEMPL>
func static_data(const Big_Int_T* data) -> const T*
{
    if constexpr(scap != 0)
        return data->static_data_;
    else
        return nullptr;
}

template <BIG_INT_TEMPL>
func alloc(Big_Int_T const& list) -> Alloc const& 
{ 
    return *cast(const Alloc*) &list; 
}

template <BIG_INT_TEMPL>
func alloc(Big_Int_T* list) -> Alloc* 
{ 
    return cast(Alloc*) list; 
}

template <BIG_INT_TEMPL>
func is_invariant(const Big_Int_T& num) -> bool
{
    const bool size_inv = num.capacity >= num.size;
    const bool data_inv = (num.capacity == 0) == (num.data == nullptr);
    const bool capa_inv = num.capacity != 0 
        ? num.capacity >= num.static_capacity
        : true;

    return size_inv && capa_inv && data_inv;
}

template <BIG_INT_TEMPL>
func is_static_alloced(Big_Int_T const& num) noexcept
{
    assert(is_invariant(num));
    return num.capacity == num.static_capacity && num.has_static_storage; 
}

template <BIG_INT_TEMPL>
proc call_grow_function(const Big_Int_T& num, size_t to_fit) -> size_t
{
    return cast(size_t) Big_Int_T::grow_type::run(
        cast(size_t) to_fit, 
        cast(size_t) num.capacity, 
        cast(size_t) num.size, 
        cast(size_t) num.static_capacity,
        sizeof(T)
    );
}

namespace detail 
{
    template<BIG_INT_TEMPL>
    proc assign(Big_Int_T* to, Big_Int_T const& other) -> void
    {
        if(to->capacity < other.size)
        {
            dealloc_data(to);
            alloc_data(to, other.size);
        }

        copy_n(to->data, other.data, other.size);
        to->size = other.size;
    }

    template<BIG_INT_TEMPL>
    proc swap(Big_Int_T* left, Big_Int_T* right) noexcept -> void 
    {
        using Big_Int = Big_Int_T;
        constexpr let transfer_static = [](Big_Int* left, Big_Int* right)
        {
            for (size_t i = 0; i < left->static_capacity; i++)
                std::swap(*static_data(left), *static_data(right));
        };

        constexpr let transfer_half_static = [=](Big_Int* from, Big_Int* to)
        {
            copy_n(static_data(to), static_data(from), to->static_capacity);

            from->data = to->data;
            to->data = static_data(to);
        };

        if (is_static_alloced(*left) && is_static_alloced(*right))
            transfer_static(left, right);
        else if (is_static_alloced(*left))
            transfer_half_static(left, right);
        else if (is_static_alloced(*right))
            transfer_half_static(right, left);
        else
            std::swap(left->data, right->data);

        std::swap(left->size, right->size);
        std::swap(left->capacity, right->capacity);

        if constexpr(left->has_stateful_alloc)
            std::swap(*alloc(left), *alloc(right));
    }

    template<BIG_INT_TEMPL>
    proc alloc_data(Big_Int_T* num, size_t capacity) -> void 
    {
        if(capacity <= num->static_capacity)
        {
            if(num->capacity < num->static_capacity)
            {
                num->capacity = num->static_capacity;
                num->data = static_data(num);
            }

            return;
        }

        num->capacity = call_grow_function(*num, capacity);
        num->data = allocate<T>(alloc(num), num->capacity * sizeof(T), DEF_ALIGNMENT<T>);
    }

    template<BIG_INT_TEMPL>
    proc dealloc_data(Big_Int_T* num) -> void 
    {
        if(num->capacity > num->static_capacity)
            deallocate<T>(alloc(num), num->data, num->capacity * sizeof(T), DEF_ALIGNMENT<T>);
    }

    //can be used for arbitrary growing/shrinking of data
    // when called on uninit Big_Int acts as alloc_data
    // when called with new_capacity = 0 acts as dealloc_data
    // Destroys elements when shrinking but does not construct new ones when growing 
    //  (because doesnt know how to)
    template<BIG_INT_TEMPL>
    proc set_capacity(Big_Int_T* num, size_t new_capacity) -> void 
    {
        T* new_data = nullptr;

        if(new_capacity <= num->static_capacity)
        {
            new_data = static_data(num);
            new_capacity = num->static_capacity;
        }
        else
        {
            let align = DEF_ALIGNMENT<T>;
            let resize_res = action<T>(alloc(num), Allocator_Actions::RESIZE, num->data, num->capacity * sizeof(T), new_capacity * sizeof(T), align, align, nullptr);
            if(resize_res.action_exists && resize_res.ptr != nullptr)
            {
                num->data = resize_res.ptr;
                num->capacity = new_capacity;
                return;
            }
            else
                new_data = allocate<T>(alloc(num), new_capacity * sizeof(T), DEF_ALIGNMENT<T>);
        }

        let copy_to = std::min(num->size, new_capacity);
        copy_n(new_data, num->data, copy_to);

        dealloc_data(num);

        num->data = new_data;
        num->capacity = new_capacity;
    }

    template<BIG_INT_TEMPL, typename T_>
    proc to_owned(Big_Int_T* num, span<T_> const& other) -> void 
    {
        alloc_data(num, other.size());
        copy_n(num->data, other.data(), other.size());
        num->size = other.size();
    }
}


template <BIG_INT_TEMPL>
proc reserve(Big_Int_T* num, size_t to_fit) -> bool
{
    assert(is_invariant(*num));

    if (num->capacity >= to_fit)
        return false;

    size_t realloc_to = 0;
    if (to_fit <= num->static_capacity)
        realloc_to = num->static_capacity;
    else
        realloc_to = call_grow_function(*num, to_fit);

    assert(realloc_to >= to_fit);
    detail::set_capacity(num, realloc_to);

    assert(is_invariant(*num));
    return true;
}

template <BIG_INT_TEMPL>
proc clear(Big_Int_T* num)
{
    assert(is_invariant(*num));
    num.size = 0;
}

template <BIG_INT_TEMPL>
func empty(Big_Int_T const& num) noexcept
{
    assert(is_invariant(num));
    return num.size == 0;
}

template <BIG_INT_TEMPL>
func is_empty(Big_Int_T const& num) noexcept
{
    assert(is_invariant(num));
    return num.size == 0;
}

template <BIG_INT_TEMPL>
proc push(Big_Int_T* num, no_infer(T) what) -> T*
{
    assert(is_invariant(*num));

    reserve(num, num->size + 1);
    num->data[num->size] = what;
    num->size++;

    assert(is_invariant(*num));
    return num->data + num->size - 1;
}

template <BIG_INT_TEMPL>
proc pop(Big_Int_T* num) -> T
{
    assert(is_invariant(*num));
    assert(num->size != 0);

    num->size--;

    return num->data[num->size];
}

template <BIG_INT_TEMPL>
proc resize(Big_Int_T* num, size_t to) -> void
{
    assert(is_invariant(*num));
    assert(0 <= to);

    let prev_size = num->size;
    reserve(num, to);
    num->size = to;

    if(prev_size < to)
        null_n(num->data + prev_size, to - prev_size);

    assert(is_invariant(*num));
}

template <BIG_INT_TEMPL>
proc resize(Big_Int_T* num, size_t to, no_infer(T) fill_with) -> void
{
    assert(is_invariant(*num));
    assert(0 <= to);

    reserve(num, to);
    num->size = to;

    for (size_t i = num->size; i < to; i++)
        (*num)[i] = fill_with;

    assert(is_invariant(*num));
}

