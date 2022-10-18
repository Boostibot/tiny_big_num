// arbitrary_precision.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <concepts>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <span>

#define let const auto
#define mut auto
#define pure [[nodiscard]]
#define proc constexpr auto
#define func pure constexpr auto
#define runtime_proc pure auto
#define runtime_func auto
#define discard cast(void)

#define cast(...) (__VA_ARGS__)
#define transmute(...) *cast(__VA_ARGS__*) cast(void*) &
#define maybe_unused [[maybe_unused]] 		
#define address_alias [[no_unique_address]]

using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

using i8 = std::int8_t;
using i16 = std::int16_t;
using i32 = std::int32_t;
using i64 = std::int64_t;

using b8 = bool;
using b16 = std::uint16_t;
using b32 = std::uint32_t;
using b64 = std::uint64_t;

using f32 = float;
using f64 = double;

using byte = u8;
using std::size_t;
using cstring = const char*;

using std::integral;
using std::span;
using std::move;

template<class T>
concept non_void = (std::is_same_v<T, void> == false);


template<typename T>
func max(T first) -> T
{   
    return first;
}

template<typename T>
func min(T first) -> T
{   
    return first;
}

template<typename T, typename ...Ts>
func max(T first, Ts... values)
    requires (std::convertible_to<Ts, T> && ...)
{
    let rest_max = max(values...);
    if(rest_max > first)
        return cast(T) rest_max;
    else
        return cast(T) first;
}

template<typename T, typename ...Ts>
func min(T first, Ts... values...)
    requires (std::convertible_to<Ts, T> && ...)
{
    let rest_max = max(values...);
    if(rest_max < first)
        return cast(T) rest_max;
    else
        return cast(T) first;
}

template<typename T_>
struct Id
{
    using Type = T_;
};

#define no_infer(...) typename Id<__VA_ARGS__>::Type

template<integral T>
func div_round_up(T value, no_infer(T) to_multiple_of) -> auto
{
    return (value + to_multiple_of - 1) / to_multiple_of;
}



namespace Allocator_Actions
{
    template <typename T>
    struct Result
    {
        bool action_exists = false;
        T* ptr = nullptr;
    };

    enum class Action : u32 {};

    constexpr Action DEALLOC_ALL = cast(Action) 1;
    constexpr Action RESIZE = cast(Action) 2;
}

template <typename T>
static constexpr size_t DEF_ALIGNMENT = max(
    alignof(std::max_align_t), 
    alignof(std::conditional_t<std::is_same_v<T, void>, byte, T>)
);

//STD allocator
template <typename Alloc>   
concept std_allocator = requires(Alloc alloc, size_t size)
{
    { alloc.allocate(size) } -> std::convertible_to<void*>;
    alloc.deallocate(nullptr, size);

    typename Alloc::value_type;
};

template <typename To, typename From>
func maybe_unsafe_ptr_cast(From* from)
{
    if constexpr (std::convertible_to<From*, To*>)
        return cast(To*) from;
    else
        return cast(To*) cast(void*) from;
}

template <typename T, std_allocator Alloc>
proc allocate(Alloc* alloc, size_t size, size_t align) -> T* 
{
    using value_type = typename Alloc::value_type;
    let recomputed_size = div_round_up(size * sizeof(T), sizeof(value_type));
    return maybe_unsafe_ptr_cast<T>(alloc->allocate(size));
}

template <typename T, std_allocator Alloc>
proc deallocate(Alloc* alloc, T* ptr, size_t size, size_t align) -> void 
{
    using value_type = typename Alloc::value_type;
    let recomputed_size = div_round_up(size * sizeof(T), sizeof(value_type));
    return alloc->deallocate(maybe_unsafe_ptr_cast<value_type>(ptr), recomputed_size);
}

template <typename T, std_allocator Alloc>
proc action(Alloc* alloc, 
    Allocator_Actions::Action action_type, 
    void* old_ptr, 
    size_t old_size, size_t new_size, 
    size_t old_align, size_t new_align, 
    void* custom_data = nullptr) -> Allocator_Actions::Result<T>
{
    return Allocator_Actions::Result<T>{false, nullptr};
}

template <typename Resource>
concept allocator = requires(Resource res, 
    Allocator_Actions::Action action_type, 
    void* old_ptr, 
    size_t old_size, size_t new_size, 
    size_t old_align, size_t new_align, 
    void* custom_data)
{
    //true; //@NOTE: MSVC is buggy and syntax highlighting breaks so we dont have all checks on per default
    { allocate<int>(&res, new_size, new_align) } -> std::convertible_to<int*>;
    deallocate<void>(&res, old_ptr, old_size, old_align);
    //{ action<void>(&res, action_type, old_ptr, old_size, new_size, old_align, new_align, custom_data) } -> std::convertible_to<Allocator_Actions::Result<void>>;
};

template<typename Container>
concept direct_container = requires(Container container)
{
    { container.data } -> std::convertible_to<void*>;
    { container.size } -> std::convertible_to<size_t>;
    requires(!std::is_same_v<decltype(container.data), void*>);
};

namespace std 
{
    func begin(direct_container auto& arr) noexcept {return arr.data;}
    func begin(const direct_container auto& arr) noexcept {return arr.data;}

    func end(direct_container auto& arr) noexcept {return arr.data + arr.size;}
    func end(const direct_container auto& arr) noexcept {return arr.data + arr.size;}

    func cbegin(const direct_container auto& arr) noexcept {return arr.data;}
    func cend(const direct_container auto& arr) noexcept {return arr.data + arr.size;}

    func size(const direct_container auto& arr) noexcept {return arr.size;}
}

template <integral T>
proc copy_n(T* to, const T* from, size_t count) -> void
{
    if (std::is_constant_evaluated())
    {
        for(size_t i = 0; i < count; i++)
            to[i] = from[i];
    }
    else
        memcpy(to, from, count * sizeof(T));
}

template <integral T>
proc null_n(T* to, const T* from, size_t count) -> void
{
    if (std::is_constant_evaluated())
    {
        for(size_t i = 0; i < count; i++)
            to[i] = 0;
    }
    else
        memset(to, 0, count * sizeof(T));
}

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

    template<integral T, size_t scap>
    func static_data(Big_Int_Data<T, scap>* data) -> T*
    {
        if constexpr(scap != 0)
            return data->static_data_;
        else
            return nullptr;
    }

    template<integral T, size_t scap>
    func static_data(const Big_Int_Data<T, scap>* data) -> const T*
    {
        if constexpr(scap != 0)
            return data->static_data_;
        else
            return nullptr;
    }
}


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

template <
    integral T, 
    size_t static_capacity_, 
    std_allocator Alloc_, 
    grower Grow
>
    struct Big_Int_;

#define Big_Int_TEMPL class T, size_t scap, class Alloc, class Grow
#define Big_Int_T Big_Int_<T, scap, Alloc, Grow>

namespace detail 
{
    template<Big_Int_TEMPL>
    proc assign(Big_Int_T* to, Big_Int_T const& other) -> void;

    template<Big_Int_TEMPL>
    proc swap(Big_Int_T* left, Big_Int_T* right) noexcept -> void;

    template<Big_Int_TEMPL>
    proc alloc_data(Big_Int_T* num, size_t capacity) -> void;

    template<Big_Int_TEMPL>
    proc dealloc_data(Big_Int_T* num) -> void;

    template<Big_Int_TEMPL, typename T_>
    proc to_owned(Big_Int_T* num, span<T_> const& other) -> void;
}

template <
    integral T_, 
    size_t static_capacity_ = 2,
    std_allocator Alloc_ = std::allocator<T_>, 
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

    constexpr Big_Int_(T val, Alloc alloc = Alloc()) 
        : Alloc{move(alloc)} 
    { 
        detail::alloc_data(this, 1);
        this->data[0] = val;
        this->size = 1;
    }

    explicit constexpr Big_Int_(slice_type slice, Alloc alloc = Alloc()) 
        : Alloc{move(alloc)} { detail::to_owned(this, slice); }

    explicit constexpr Big_Int_(const_slice_type slice, Alloc alloc = Alloc()) 
        : Alloc{move(alloc)} { detail::to_owned(this, slice); }

    constexpr Big_Int_(T* data, size_t size, size_t capacity, Alloc alloc, const T (&static_data)[static_capacity]) noexcept 
        requires (has_static_storage)
        : Alloc{move(alloc)}, Big_Int_Data{data, size, capacity} 
    {
        copy_n(static_data(this), static_data. static_capacity);
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
        { return slice_type{this->data, this->size}; }

    constexpr operator const_slice_type() const noexcept 
        { return const_slice_type{this->data, this->size}; }

};

namespace std 
{
    template <Big_Int_TEMPL>
    proc swap(Big_Int_T& num1, Big_Int_T& num2)
    {
        return detail::swap(&num1, &num2);
    }
}

template <Big_Int_TEMPL>
func alloc(Big_Int_T const& list) -> Alloc const& { return *cast(const Alloc*) &list; }

template <Big_Int_TEMPL>
func alloc(Big_Int_T* list) -> Alloc* { return cast(Alloc*) list; }

template <Big_Int_TEMPL>
func is_invariant(const Big_Int_T& num) -> bool
{
    const bool size_inv = num.capacity >= num.size;
    const bool data_inv = (num.capacity == 0) == (num.data == nullptr);
    const bool capa_inv = num.capacity != 0 
        ? num.capacity >= num.static_capacity
        : true;

    return size_inv && capa_inv && data_inv;
}

template <Big_Int_TEMPL>
func is_static_alloced(Big_Int_T const& num) noexcept
{
    assert(is_invariant(num));
    return num.capacity == num.static_capacity && num.has_static_storage; 
}

template <Big_Int_TEMPL>
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
    template<Big_Int_TEMPL>
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

    template<Big_Int_TEMPL>
    proc swap(Big_Int_T* left, Big_Int_T* right) noexcept -> void 
    {
        using Big_Int = Big_Int_T;
        constexpr let transfer_static = [](Big_Int* left, Big_Int* right)
        {
            size_t to_size = min(left->size, right->size);
            for (size_t i = 0; i < to_size; i++)
                std::swap((*left)[i], (*right)[i]);
        };

        constexpr let transfer_half_static = [=](Big_Int* from, Big_Int* to)
        {
            transfer_static(from, to, 0, from->size);

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

    template<Big_Int_TEMPL>
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

    template<Big_Int_TEMPL>
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
    template<Big_Int_TEMPL>
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

    template<Big_Int_TEMPL, typename T_>
    proc to_owned(Big_Int_T* num, span<T_> const& other) -> void 
    {
        alloc_data(num, other.size());
        copy_n(num->data, other.data(), other.size());
        num->size = other.size();
    }
}


template <Big_Int_TEMPL>
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

template <Big_Int_TEMPL>
proc clear(Big_Int_T* num)
{
    assert(is_invariant(*num));
    num.size = 0;
}

template <Big_Int_TEMPL>
func empty(Big_Int_T const& num) noexcept
{
    assert(is_invariant(num));
    return num.size == 0;
}

template <Big_Int_TEMPL>
func is_empty(Big_Int_T const& num) noexcept
{
    assert(is_invariant(num));
    return num.size == 0;
}

template <Big_Int_TEMPL>
func operator +=(Big_Int_T& num, Big_Int_T const& num2)
{

}


int main()
{
    Big_Int_<u64> num;
    reserve(&num, 20);

    let x = div_round_up(5, 5.0);

    u64 arr[] = {1, 2, 3};
    span s{arr, 3};
    span<const u64> sc{arr, 3};

    { Big_Int_<u64> num(s); }
    { Big_Int_<u64> num(sc); }
    { 
        Big_Int_<u64> num = 10; 
        for(auto val: num)    
            std::cout << val;
    }
    { Big_Int_<u64> num = Big_Int_<u64>(std::allocator<u64>()); }

    std::cout << "\nOK\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
