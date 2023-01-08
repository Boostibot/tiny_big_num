#include <cstdint>

using u64 = uint64_t;
using u128 = __uint128_t;

// proof of concept, actual size is size*64bits
template<int size> struct Batch
{
    u64 data[size];
};

//clang
inline u64 add_with_carry(u64 &a, const u64 &b, const u64 &c)
{
    u128 s = u128(a) + b + c;
    a = (u64) s;
    return s >> 64;
}

//gnu
inline u64 add_with_carry(u64 carry, u64 a, u64 b, u64* out)
{
    u64 res = a + b + carry;
    *out = res;
    return (u64) (res < a);
}

//msvc c = _addcarry_u64(c, left0, right0, out + 0);
void test(Batch<3> &a, const Batch<3> &b_)
{
    const Batch<3> b = b_; // need to copy b_ into local temporary
    u64 c0 = add_with_carry(a.data[0], b.data[0], 0);
    u64 c1 = add_with_carry(a.data[1], b.data[1], c0);
    u64 c2 = add_with_carry(a.data[2], b.data[2], c1);
}


//gcc & msvc -> only one add in place can be the same
//clang -> two adds

//clang & gcc -> only one mul for in place and out of place
//msvc -> two muls (just the inline 10 lines suffices

// generic solution (no border case handling for size=0):

carry_in = mul_carry(carry_in, out->data[0], b, &out->data[0]);
carry_in = mul_carry(carry_in, out->data[1], b, &out->data[1]);
carry_in = mul_carry(carry_in, out->data[2], b, &out->data[2]);
carry_in = mul_carry(carry_in, out->data[3], b, &out->data[3]);
carry_in = mul_carry(carry_in, out->data[4], b, &out->data[4]);
carry_in = mul_carry(carry_in, out->data[5], b, &out->data[5]);
carry_in = mul_carry(carry_in, out->data[6], b, &out->data[6]);
carry_in = mul_carry(carry_in, out->data[7], b, &out->data[7]);
carry_in = mul_carry(carry_in, out->data[8], b, &out->data[8]);
carry_in = mul_carry(carry_in, out->data[9], b, &out->data[9]);
template<int size>
u64 addTo2(Batch<size> &a, const Batch<size> b, u64 carry_in)
{
    u128 c = u128(a.data[0]) + b.data[0] + carry_in;
    a.data[0] = c;

    for(int i=1; i<size; ++i)
    {
        c = u128(a.data[i]) + b.data[i] + (c >> 64);
        a.data[i] = c;
    }

    return (c >> 64);
}

template<int size>
void add(Batch<size>* out, const Batch<size> a, const Batch<size> b, u64 carry_in)
{
    u128 c = u128(a.data[0]) + b.data[0] + carry_in;
    for(int i=1; i<size; ++i)
    {
        c = u128(a.data[i]) + b.data[i] + (c >> 64);
        out->data[i] = c;
    }

    return (c >> 64);
}


void test10_2(Batch<10> &a, const Batch<10>& b)
{
    addTo2(a, b, 0);
}
template<int size, int i>
u64 add_gcc_(Batch<size> &a, const Batch<size>& b)
{
    a.data[i] += b.data[i];
    a.data[i] += a.data[i - 1] < b.data[i - 1];
    if constexpr(i != size - 1)
        return add_gcc_<size, i + 1>(a, b);
    else
    {
        return (u64) (a.data[i] < b.data[i]);
    }
}
template<int size>
u64 add_gcc(Batch<size> &a, const Batch<size>& b, u64 carry)
{
    a.data[0] += b.data[0] + carry;
    return add_gcc_<size, 1>(a, b);
}


u64 add_gcc_8(Batch<8>& a, Batch<8> b, u64 carry)
{
    return add_gcc<8>(a, b, carry);
}
u64 manual_carry8_out_in(
    u64* out, u64 carry,
    u64 left0, u64 left1, u64 left2, u64 left3, u64 left4, u64 left5, u64 left6, u64 left7, 
    u64 right0, u64 right1, u64 right2, u64 right3, u64 right4, u64 right5, u64 right6, u64 right7)
{
    left0 += right0 + carry;
    left1 += right1;
    left1 += left0 < right0;
    left2 += right2;
    left2 += left1 < right1;
    left3 += right3;
    left3 += left2 < right2;

    left4 += right4;
    left4 += left3 < right3;
    left5 += right5;
    left5 += left4 < right4;
    left6 += right6;
    left6 += left5 < right5;
    left7 += right7;
    left7 += left6 < right6;

    out[0] = left0;
    out[1] = left1;
    out[2] = left2;
    out[3] = left3;
    out[4] = left4;
    out[5] = left5;
    out[6] = left6;
    out[7] = left7;

    return (u64) (left7 < right7);
}

#if 0
template<int size, int i = 0>
u64 mul(Batch<size>* out, const Batch<size> a, u64 b, u64 carry_in)
{
    if constexpr (i == size)
        return carry_in;
    else
    {
        u128 c = u128(a.data[i]) * b;
        u128 over = (c >> 64) + carry_in;
        out->data[i] = (over >> 64);
        carry_in = over;

        return mul<size, i + 1>(out, a, b, carry_in);
    }
}
#else
template<int size>
void mul(Batch<size>* out, const Batch<size> a, u64 b, u64 carry_in)
{
    u128 c = u128(a.data[0]) * b;
    u128 over = (c >> 64) + carry_in;
    out->data[0] = (over >> 64);
    carry_in = over;

    for(int i=1; i<size; ++i)
    {
        c = u128(a.data[i]) * b;
        over = (c >> 64) + carry_in;
        out->data[i] = (over >> 64);
        carry_in = over;
    }
}
#endif


u64 add_10(Batch<10>* out, const Batch<10> left, const Batch<10> right, u64 carry_in)
{
    carry_in = add_carry(carry_in, left.data[0], right.data[0], out->data + 0);
    carry_in = add_carry(carry_in, left.data[1], right.data[1], out->data + 1);
    carry_in = add_carry(carry_in, left.data[2], right.data[2], out->data + 2);
    carry_in = add_carry(carry_in, left.data[3], right.data[3], out->data + 3);
    carry_in = add_carry(carry_in, left.data[4], right.data[4], out->data + 4);
    carry_in = add_carry(carry_in, left.data[5], right.data[5], out->data + 5);
    carry_in = add_carry(carry_in, left.data[6], right.data[6], out->data + 6);
    carry_in = add_carry(carry_in, left.data[7], right.data[7], out->data + 7);
    return carry_in;
}

Batch<2> mul_overflow(u64 l, u64 r)
{
    u128 res = l;
    res *= r;
    Batch<2> out;
    out.data[0] = (u64) res;
    out.data[1] = res >> 64;
    return out;
}


void mul_10(Batch<10>* out, const Batch<10> a, u64 b, u64 carry_in)
{
    mul(out, a, b, carry_in);
}

#include <iostream>

int main()
{
    u64 x = -1;
    auto res = mul_overflow(x, 4);
    x *= 7;
    std::cout << res.data[1] << std::endl;


}

#include <stdint.h>
#include <inttypes.h>

using u64 = uint64_t;

#ifdef _MSC_VER
#include <intrin.h>

void manual_carry8_out_intrin(
    u64* out,
    u64 left0, u64 left1, u64 left2, u64 left3, u64 left4, u64 left5, u64 left6, u64 left7, 
    u64 right0, u64 right1, u64 right2, u64 right3, u64 right4, u64 right5, u64 right6, u64 right7)
{
    unsigned char c = 0;
    c = _addcarry_u64(c, left0, right0, out + 0);
    c = _addcarry_u64(c, left1, right1, out + 1);
    c = _addcarry_u64(c, left2, right2, out + 2);
    c = _addcarry_u64(c, left3, right3, out + 3);
    c = _addcarry_u64(c, left4, right4, out + 4);
    c = _addcarry_u64(c, left5, right5, out + 5);
    c = _addcarry_u64(c, left6, right6, out + 6);
    c = _addcarry_u64(c, left7, right7, out + 7);
}

template<int size> struct Batch
{
    u64 data[size];
};

template<int size>
struct Counter{};

template<int i, int size>
inline u64 mul_intrin(Batch<size>* out, const Batch<size> a, u64 b, u64 carry_in)
{
    u64 r2 = 0;
    u64 r1 = _umul128(a.data[i], b, &r2);

    unsigned char c = 0;
    c = _addcarry_u64(c, (u64) r1, carry_in, out->data + i);
    carry_in = c + (u64) r2;

    return carry_in;
}


template<int size>
inline u64 mul_intrin(Batch<size>* out, const Batch<size> a, u64 b, u64 carry_in, Counter<size>)
{
    return carry_in;
}
#endif

u64 mul(Batch<10>* out, const Batch<10> a, u64 b, u64 carry_in)
{
    carry_in = mul_intrin<0>(out, a, b, carry_in);
    carry_in = mul_intrin<1>(out, a, b, carry_in);
    carry_in = mul_intrin<2>(out, a, b, carry_in);
    carry_in = mul_intrin<3>(out, a, b, carry_in);
    carry_in = mul_intrin<4>(out, a, b, carry_in);
    carry_in = mul_intrin<5>(out, a, b, carry_in);
    carry_in = mul_intrin<6>(out, a, b, carry_in);
    carry_in = mul_intrin<7>(out, a, b, carry_in);
    carry_in = mul_intrin<8>(out, a, b, carry_in);
    carry_in = mul_intrin<9>(out, a, b, carry_in);

    return carry_in;
}


//NEW
#include <stdint.h>
#include <inttypes.h>

using u64 = uint64_t;

template<int size> struct Batch
{
    u64 data[size];
};

#if defined(__clang__)
using u128 = __uint128_t;

inline u64 add_carry(u64 carry, u64 a, u64 b, u64* out)
{
    u128 s = u128(a) + b + carry;
    *out = (u64) s;
    return s >> 64;
}

template<int size>
u64 add(Batch<size>* out, const Batch<size> a, const Batch<size> b, u64 carry_in)
{
    u128 c = u128(a.data[0]) + b.data[0] + carry_in;
    for(int i=1; i<size; ++i)
    {
        c = u128(a.data[i]) + b.data[i] + (c >> 64);
        out->data[i] = c;
    }

    return (c >> 64);
}

void add_10_2_in_place(Batch<10> &a, const Batch<10>& b)
{
    add(&a, a, b, 0);
}
#elif defined(__GNUC__)
using u128 = __uint128_t;

inline u64 add_carry(u64 carry, u64 a, u64 b, u64* out)
{
    u64 res = a + b + carry;
    *out = res;
    return (u64) (res < a);
}
#endif

#ifdef __GNUC__
u64 mul_carry(u64 carry, u64 a, u64 b, u64* out)
{
    u128 over = u128(a) * b + carry;
    *out = (u64) over;
    carry = (over >> 64);
    return over;
}
#endif

#ifdef _MSC_VER
#include <intrin.h>

inline u64 add_carry(u64 carry, u64 a, u64 b, u64* out)
{
    return (u64) _addcarry_u64((unsigned char) carry, a, b, out);
}

inline u64 mul_carry(u64 carry, u64 a, u64 b, u64* out)
{
    u64 r2 = 0;
    u64 r1 = _umul128(a, b, &r2);

    unsigned char c = 0;
    c = _addcarry_u64(c, r1, carry, out);
    carry = c + (u64) r2;
    return carry;
}

u64 mul10_in_place_msvc(Batch<10>* out, u64 b, u64 carry_in)
{
    carry_in = mul_carry(carry_in, out->data[0], b, &out->data[0]);
    carry_in = mul_carry(carry_in, out->data[1], b, &out->data[1]);
    carry_in = mul_carry(carry_in, out->data[2], b, &out->data[2]);
    carry_in = mul_carry(carry_in, out->data[3], b, &out->data[3]);
    carry_in = mul_carry(carry_in, out->data[4], b, &out->data[4]);
    carry_in = mul_carry(carry_in, out->data[5], b, &out->data[5]);
    carry_in = mul_carry(carry_in, out->data[6], b, &out->data[6]);
    carry_in = mul_carry(carry_in, out->data[7], b, &out->data[7]);
    carry_in = mul_carry(carry_in, out->data[8], b, &out->data[8]);
    carry_in = mul_carry(carry_in, out->data[9], b, &out->data[9]);
    return carry_in;
}

#endif

u64 add_10(Batch<10>* out, const Batch<10> left, const Batch<10> right, u64 carry_in)
{
    carry_in = add_carry(carry_in, left.data[0], right.data[0], out->data + 0);
    carry_in = add_carry(carry_in, left.data[1], right.data[1], out->data + 1);
    carry_in = add_carry(carry_in, left.data[2], right.data[2], out->data + 2);
    carry_in = add_carry(carry_in, left.data[3], right.data[3], out->data + 3);
    carry_in = add_carry(carry_in, left.data[4], right.data[4], out->data + 4);
    carry_in = add_carry(carry_in, left.data[5], right.data[5], out->data + 5);
    carry_in = add_carry(carry_in, left.data[6], right.data[6], out->data + 6);
    carry_in = add_carry(carry_in, left.data[7], right.data[7], out->data + 7);
    return carry_in;
}




u64 add_10_in_place(Batch<10>* out, const Batch<10> right, u64 carry_in)
{
    // for(int i=0; i<10; ++i)
    // {
    //     u64 a = out->data[i];
    //     u64 b = right.data[i];
    //     carry_in = add_carry(carry_in, a, b, &a);
    //     out->data[i] = a;
    // }
    // return carry_in;

    carry_in = add_carry(carry_in, out->data[0], right.data[0], out->data + 0);
    carry_in = add_carry(carry_in, out->data[1], right.data[1], out->data + 1);
    carry_in = add_carry(carry_in, out->data[2], right.data[2], out->data + 2);
    carry_in = add_carry(carry_in, out->data[3], right.data[3], out->data + 3);
    carry_in = add_carry(carry_in, out->data[4], right.data[4], out->data + 4);
    carry_in = add_carry(carry_in, out->data[5], right.data[5], out->data + 5);
    carry_in = add_carry(carry_in, out->data[6], right.data[6], out->data + 6);
    carry_in = add_carry(carry_in, out->data[7], right.data[7], out->data + 7);
    return carry_in;
}

inline __forceinline u64 mul10(Batch<10>* out, const Batch<10> a, u64 b, u64 carry_in)
{
    carry_in = mul_carry(carry_in, a.data[0], b, &out->data[0]);
    carry_in = mul_carry(carry_in, a.data[1], b, &out->data[1]);
    carry_in = mul_carry(carry_in, a.data[2], b, &out->data[2]);
    carry_in = mul_carry(carry_in, a.data[3], b, &out->data[3]);
    carry_in = mul_carry(carry_in, a.data[4], b, &out->data[4]);
    carry_in = mul_carry(carry_in, a.data[5], b, &out->data[5]);
    carry_in = mul_carry(carry_in, a.data[6], b, &out->data[6]);
    carry_in = mul_carry(carry_in, a.data[7], b, &out->data[7]);
    carry_in = mul_carry(carry_in, a.data[8], b, &out->data[8]);
    carry_in = mul_carry(carry_in, a.data[9], b, &out->data[9]);

    return carry_in;
}

u64 mul10_in_place(Batch<10>* out, u64 b, u64 carry_in)
{
    return mul10(out, *out, b, carry_in);
}