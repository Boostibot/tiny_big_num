#pragma once

#include <concepts>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <span>
#include <ranges>

#define let const auto
#define mut auto
#define pure [[nodiscard]]
#define proc constexpr auto
#define func pure constexpr auto
#define runtime_proc pure auto
#define runtime_func auto

#define cast(...) (__VA_ARGS__)
#define transmute(...) *cast(__VA_ARGS__*) cast(void*) &
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

template<class T1, class T2>
concept same = std::is_same_v<T1, T2>;

template<class T>
concept non_void = (!same<T, void>);


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
func max(T first, Ts... values) -> T
    requires (std::convertible_to<Ts, T> && ...)
{
    let rest_max = max(values...);
    if(rest_max > first)
        return cast(T) rest_max;
    else
        return cast(T) first;
}

template<typename T, typename ...Ts>
func min(T first, Ts... values...) -> T
    requires (std::convertible_to<Ts, T> && ...)
{
    let rest_max = max(values...);
    if(rest_max < first)
        return cast(T) rest_max;
    else
        return cast(T) first;
}

template<typename T>
using No_Infer = std::type_identity<T>::type;
#define no_infer(...) No_Infer<__VA_ARGS__> 

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
    alignof(std::conditional_t<non_void<T>, T, byte>)
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
concept allocator = true; 

//template <typename Resource>
//concept allocator = requires(Resource res, 
//    Allocator_Actions::Action action_type, 
//    void* old_ptr, 
//    size_t old_size, size_t new_size, 
//    size_t old_align, size_t new_align, 
//    void* custom_data)
//{
//    //true; //@NOTE: MSVC is buggy and syntax highlighting breaks so we dont have all checks on per default
//    //{ allocate<int>(&res, new_size, new_align) } -> std::convertible_to<int*>;
//    //deallocate<void>(&res, old_ptr, old_size, old_align);
//    //{ action<void>(&res, action_type, old_ptr, old_size, new_size, old_align, new_align, custom_data) } -> std::convertible_to<Allocator_Actions::Result<void>>;
//};

template<typename Container>
concept direct_container = requires(Container container)
{
    { container.data } -> std::convertible_to<void*>;
    { container.size } -> std::convertible_to<size_t>;
    requires(!same<decltype(container.data), void*>);
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
proc null_n(T* to, size_t count) -> void
{
    if (std::is_constant_evaluated())
    {
        for(size_t i = 0; i < count; i++)
            to[i] = 0;
    }
    else
        memset(to, 0, count * sizeof(T));
}


template <typename T, size_t size_>
struct Array
{
    static constexpr size_t size = size_;
    T data[size > 0 ? size : 1];

    func& operator[](size_t index) const noexcept { assert(index < size && "index out of range"); return data[index]; }
    func& operator[](size_t index) noexcept       { assert(index < size && "index out of range"); return data[index]; }
};