#include <iostream>
#include <vector>
#include <concepts>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <span>
#include <limits>

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
    //requires (std::convertible_to<Ts, T> && ...)
{
    let rest_max = max(values...);
    if(rest_max > first)
        return cast(T) rest_max;
    else
        return cast(T) first;
}

template<typename T, typename ...Ts>
func min(T first, Ts... values...)
    //requires (std::convertible_to<Ts, T> && ...)
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
    std_allocator Alloc_, 
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

    reserve(num, to);
    num->size = to;

    if(num->size < to)
        null_n(num->data + num->size, to - num->size);

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



template <typename T>
constexpr size_t BIT_SIZE = sizeof(T) * CHAR_BIT;
template <typename T>
constexpr size_t HALF_BIT_SIZE = BIT_SIZE<T> / 2;

template <integral T>
constexpr T LOW_MASK = cast(T)(cast(T)(-1) >> HALF_BIT_SIZE<T>);
template <integral T>
constexpr T HIGH_MASK = cast(T)(cast(T)(-1) << HALF_BIT_SIZE<T>);
template <integral T>
constexpr T HIGH_1 = cast(T)(cast(T)(1) << HALF_BIT_SIZE<T>); 

template <integral T>
func high_bits(T value) -> T {
    return value >> HALF_BIT_SIZE<T>;
};  

template <integral T>
func low_bits(T value) -> T {
    return value & LOW_MASK<T>;
};  

template <integral T>
func combine_bits(T low, T high) -> T {
    return (low & LOW_MASK<T>) | (high << HALF_BIT_SIZE<T>);
};

namespace detail
{
    enum class Additive_Mode 
    {
        ADD,
        SUB
    };

    template <Additive_Mode mode, integral T>
    proc add_or_sub_assign_carry(span<T>* to, span<const T> from) -> T
    {
        assert(to->size() >= from.size());

        let add_sub_half = [](T last, T left, T right){
            if constexpr(mode == Additive_Mode::ADD)
                return left + right + high_bits(last);
            else
                return left - right + high_bits(last);
        };  

        T last = 0;
        for(size_t i = 0; i < from.size(); i++)
        {
            const T from_val = from[i];
            const T to_val = (*to)[i];

            const T res_low = add_sub_half(last, low_bits(from_val), low_bits(to_val));
            const T res_high = add_sub_half(res_low, high_bits(from_val), high_bits(to_val));

            const T res = combine_bits(res_low, res_high);

            (*to)[i] = res;
            last = res;
        }

        return high_bits(last) > 0 ? 1 : 0;
    }
}

template <integral T>
func add_assign_carry(span<T>* to, span<const T> from) -> T
{
    assert(to->size() >= from.size());
    
    T carry = detail::add_or_sub_assign_carry<detail::Additive_Mode::ADD>(to, from);
    if(carry == 0)
        return carry;

    //iterate untill carry is exhausted 
    for (size_t i = from.size(); i < to->size(); i++)
    {
        (*to)[i] += 1;

        if((*to)[i] != 0)
        {
            carry = 0;
            break;
        }
    }

    return carry;
}


template <integral T>
func sub_assign_carry(span<T>* to, span<const T> from) -> T
{
    assert(to->size() >= from.size());

    T carry = detail::add_or_sub_assign_carry<detail::Additive_Mode::SUB>(to, from);
    if(carry == 0)
        return carry;

    //iterate untill carry is exhausted 
    for (size_t i = from.size(); i < to->size(); i++)
    {
        //if wont underflow
        if((*to)[i] != 0)
        {
            (*to)[i] -= 1;
            carry = 0;
            break;
        }

        (*to)[i] -= 1;
    }

    return carry;
}

template <integral T>
func zeros_from_index(span<T> num) -> span<T>
{
    size_t to = num.size();
    for(; to-- > 0;)
        if(num[to] != 0)
            break;

    return to + 1;
}

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <chrono>

template <BIG_INT_TEMPL>
runtime_proc print_digits(Big_Int_T const& num)
{
    for(let word : num)
        std::cout << cast(u64) word;
}


template <typename Fn, typename To = void>
concept thunk = requires(Fn fn)
{
    { fn() } -> std::convertible_to<To>;
};




#if defined(__GNUC__)
    #define PERF_UNUSED __attribute__((unused))
    #define ALWAYS_INLINE __attribute__((always_inline))
    #define inline_assembly asm
#elif defined(_MSC_VER) && !defined(__clang__)
    #define PERF_UNUSED
    #define ALWAYS_INLINE __forceinline
    #define inline_assembly __asm
#else
    #define PERF_UNUSED
    #define ALWAYS_INLINE
    #define inline_assembly asm
#endif


inline ALWAYS_INLINE proc no_optim(void* value) -> void 
{
    if(!std::is_constant_evaluated())
    {
        // Clang doesn't like the 'X' constraint on `value` and certain GCC versions
        // don't like the 'g' constraint. Attempt to placate them both.
        /*#if defined(__clang__)
        inline_assembly volatile("" : : "g"(value) : "memory");
        #else
        inline_assembly volatile("" : : "i,r,m"(value) : "memory");
        #endif*/
    }
}
// Force the compiler to flush pending writes to global memory. Acts as an
// effective read/write barrier
inline ALWAYS_INLINE proc flush_writes() -> void {

    /*if(!std::is_constant_evaluated())
    {
        inline_assembly volatile("" : : : "memory");
    }*/
}



template <thunk Fn>
func ellapsed_time(Fn fn) -> u64
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    let t1 = high_resolution_clock::now();
    fn();
    let t2 = high_resolution_clock::now();
    return duration_cast<milliseconds>(t2 - t1).count();
}



template <typename Fn>
func benchmark_iters(size_t ms, size_t warm_up_ms, Fn fn, size_t block_size = 1) -> size_t
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    let clock = [](){ 
        return high_resolution_clock::now(); 
    };
    let diff = [](auto t1, auto t2){ 
        return duration_cast<milliseconds>(t2 - t1).count(); 
    };

    size_t total_ms = ms + warm_up_ms;
    size_t iters = 0;
    let start = clock();

    while(true)
    {
        for(size_t i = 0; i < block_size; i++)
            fn();

        let elapsed = cast(size_t) diff(start, clock());
        if(elapsed > warm_up_ms)
            iters += block_size;
        if(elapsed > total_ms)
            break;
    }

    return iters;
}

struct benchmark_result
{
    size_t iters;
    size_t time;
};

template <typename Fn>
func benchmark_iters2(size_t min_time, size_t max_time, Fn fn) -> double
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    let clock = [](){ 
        return high_resolution_clock::now(); 
    };
    let diff = [](auto t1, auto t2){ 
        return duration_cast<milliseconds>(t2 - t1).count(); 
    };

    enum class Stage
    {
        Initial,
        Blocked
    };

    constexpr double acceptable_delta_ratio = 0.01;
    constexpr double perfect_delta_ratio = 0.01;
    constexpr size_t enough_samples = 50;
    constexpr size_t history_len = 8;
    double block_times[history_len + 1] = {0};
    double block_diffs_[history_len + 1] = {0.0};
    double* block_diffs = block_diffs_ + 1;

    size_t block_index = 0;

    size_t time_per_block = 50;
    //size_t total_ms = ms + warm_up_ms;
    size_t iters = 0;
    size_t samples = 0;
    size_t last_iters = 0;
    size_t block_size = 1;
    size_t to_time = min_time;
    Stage stage = Stage::Initial;


    size_t reported_iters = 0;
    double reported_time = 0;
    double max_diff = 0;

    let start = clock();
    mut last = start;
    double time_per_iter = 0.0;

    while(true)
    {
        for(size_t i = 0; i < block_size; i++)
            fn();

        iters += block_size;

        let now = clock();
        let ellapsed = cast(size_t) diff(start, now);
        if(ellapsed > to_time)
        {
            if(ellapsed > max_time)
                break;

            if(stage == Stage::Initial)
            {
                time_per_iter = cast(double) ellapsed / iters;
                //let expected_iters = optimal_time * iters / ellapsed; //== optimal_time / time_per iter
                block_size = time_per_block * iters / ellapsed; // == time_per_block / time_per_iter;
                
                stage = Stage::Blocked;
                to_time = time_per_block;
            }
            else
            {
                double current_iter_time = cast(double) diff(last, now);
                size_t current_iter_count = iters - last_iters;

                size_t prev_index = block_index % history_len;
                size_t curr_index = prev_index + 1;

                block_times[curr_index] = current_iter_time;
                double delta = block_times[curr_index] - block_times[prev_index];
                double delta2 = delta * delta;
                block_diffs[curr_index] = delta * delta;
                double max = 0;

                if(delta2 > max)
                    max = delta2;
                    

                block_index ++;
                
                if(block_index > history_len - 1)
                {
                    double max_delta2 = 0;
                    for(size_t i = 0; i < history_len - 1; i++)
                    {
                        double delta = block_times[i] - block_times[i+1];
                        double delta2 = delta * delta;
                        if(delta2 > max_delta2)
                            max_delta2 = delta2;
                    }


                    double normalized_delta = max_delta2 / (current_iter_time * current_iter_time);
                    if(normalized_delta < acceptable_delta_ratio)
                    {
                        reported_iters += current_iter_count;
                        reported_time += current_iter_time;
                        samples++;
                    }

                    if(samples > enough_samples)
                        break;
                    //if(normalized_delta < perfect_delta_ratio)
                        //break;
                }
            }

            //last = now;
            last = clock();
            last_iters = iters;
        }
    }

    if(reported_iters == 0)
    {
        reported_iters = iters;
        reported_time = cast(size_t) diff(start, clock());
    }

    std::cout << "total_time: " << diff(start, clock()) << "\n";

    return cast(double) reported_time / cast(double) reported_iters;
    //return benchmark_result{reported_iters, cast(size_t) reported_time};
}

double run_test()
{
    using Big_Int = Big_Int_<u64>;
    constexpr size_t MIN_NUM_SIZE = 10000;
    constexpr size_t MAX_NUM_SIZE = 10000;
    constexpr size_t NUM_COUNT = 20;
    constexpr size_t WARM_ITERS = 100;
    constexpr size_t ITERS = 1000'0;

    let gen_random = [](size_t min_num_size, size_t max_num_size) -> Big_Int {
        size_t size = min(max_num_size, max(min_num_size, rand()));

        Big_Int gened;
        resize(&gened, size);
        for(size_t i = 0; i < size; i++)
        {
            const u64 r = rand();
            gened[i] = r | (r >> 16) | (r >> 32) | (r >> 48);
        }

        return gened;
    };

    Big_Int to = gen_random(MAX_NUM_SIZE, MAX_NUM_SIZE);
    Big_Int nums[NUM_COUNT];
    for(size_t i = 0; i < NUM_COUNT; i++)
        nums[i] = gen_random(MIN_NUM_SIZE, MAX_NUM_SIZE);

    let lambda = [&]{
        span<u64> to_s = to;
        for(size_t i = 0; i < NUM_COUNT; i++)
        {
            span<const u64> num_s = nums[i];
            cast(void) add_assign_carry(&to_s, num_s);
        }
    };

    //std::cout << "iters:" << benchmark_iters(1000, 2000, lambda, 10) << " \n";
    //std::cout << "iters:" << benchmark_iters(1000, 1000, lambda, 10) << " \n";
    //std::cout << "iters:" << benchmark_iters(1000, 500, lambda, 10) << " \n";
    //std::cout << "iters:" << benchmark_iters(1000, 100, lambda, 10) << " \n";
    //std::cout << "iters:" << benchmark_iters(1000, 100, lambda, 10) << " \n";

    std::cout << "iters:" << benchmark_iters2(100, 1000000, lambda) << " \n";
    std::cout << "iters:" << benchmark_iters2(100, 1000000, lambda) << " \n";
    std::cout << "iters:" << benchmark_iters2(100, 1000000, lambda) << " \n";
    std::cout << "iters:" << benchmark_iters2(100, 1000000, lambda) << " \n";


    return 1.0;
}

int main()
{
    Big_Int_<u8> num1 = 0xFF;
    Big_Int_<u8> num2 = 0xA2;

    constexpr int res1 = 0xFF + 0xA2;
    constexpr int res2 = 0x1b1;


    resize(&num1, 2);
    
    span<u8> num1_s = num1;
    span<const u8> num2_s = num2;
    int a = 5;

    run_test();

    std::cout << "\nOK\n";
}