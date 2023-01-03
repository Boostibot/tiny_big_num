#pragma once

#include <stdexcept>
#include <iostream>
#include "pod_vector.h"
#include "tiny_big_num.h"

namespace big_int
{
    using namespace tiny_num;

    using imax = int64_t;
    using umax = uint64_t;
    using Digit = umax;
    using Digits = POD_Vector<Digit>;
    using Slice = tiny_num::Slice<Digit>;
    using CSlice = tiny_num::Slice<const Digit>;
    using Single_Overflow = tiny_num::Single_Overflow<Digit>;
    using Batch_Overflow = tiny_num::Batch_Overflow<Digit>;
    static constexpr Optims optims = Optims{};

    inline static Slice slice(Digits& digits) noexcept {
        return Slice{digits.data(), (size_t) digits.size()};
    }

    inline static CSlice slice(Digits const& digits) noexcept {
        return CSlice{digits.data(), (size_t) digits.size()};
    }

    inline static Slice buffer(Digits& digits) noexcept {
        return Slice{digits.data(), (size_t) digits.capacity()};
    }

    inline static CSlice buffer(Digits const& digits) noexcept {
        return CSlice{digits.data(), (size_t) digits.capacity()};
    }


    class Big_Int_Pusher
    {
        Push_Digit<umax> pusher;
        Digits digits;
        int sign = 1;

        void resize_digits(size_t size)
        {
            digits.resize(size);
            pusher.buffer = slice(digits);
        }

    public:
        Big_Int_Pusher(Big_Int_Pusher const&) = delete; 

        Big_Int_Pusher()
        {
            Slice my_buffer = slice(digits);
            pusher = push_digit_init(&my_buffer);
        }

        void set_sign(int sign) {
            this->sign = sign;
        }

        void push(umax digit, umax base)
        {
            bool ok = push_digit(&pusher, digit, base, optims);
            if(!ok)
            {
                resize_digits(pusher.num_size + 1);
                ok = push_digit(&pusher, digit, base, optims);
                assert(ok);
            }
        }

        friend class Big_Int;
    };

    class Big_Int
    {

    private:
        Digits _digits;
        Digits _aux;
        int _sign = 1;

    public:
        Big_Int() noexcept = default; 
        
        Big_Int(imax digit) 
        {
            if(digit == 0)
                return;

            if(digit < 0)
            {
                _sign = -1;
                digit = -digit;
            }

            _digits.push_back(digit);
        }

        Big_Int(Big_Int_Pusher&& pusher)
        {
            auto res = push_digit_result(&pusher.pusher, optims);
            if(!res.ok)
            {
                pusher.resize_digits(pusher.pusher.num_size + 1);
                res = push_digit_result(&pusher.pusher, optims);
                assert(res.ok);
            }

            pusher.resize_digits(res.slice.size);
            _digits = std::move(pusher.digits);
            _sign = pusher.sign;
        }

    public:
        inline Slice  slice() noexcept          { return big_int::slice(_digits); }
        inline CSlice slice() const noexcept    { return big_int::slice(_digits); }
        inline Slice  buffer() noexcept         { return big_int::buffer(_digits); }
        inline CSlice buffer() const noexcept   { return big_int::buffer(_digits); }
        inline Slice  aux() noexcept            { return big_int::buffer(_aux); }
        inline CSlice aux() const noexcept      { return big_int::buffer(_aux); }

        inline int  sign() const noexcept        { return _sign; }
        inline void set_sign(int value) noexcept { _sign = value; }
        inline void flip_sign() noexcept         { _sign *= -1; }
        
        inline size_t digit_count() const noexcept { 
            return _digits.size(); 
        }

    private:
        inline void resize_slice(size_t size) {
            _digits.resize_for_overwrite(size);
        }

        inline void grow_buffer(size_t size) {
            _digits.reserve(size);
        }

        inline void grow_aux(size_t size) {
            _aux.reserve(size);
        }

        size_t size() const noexcept {
            return _digits.size();
        }

    private:
        void add_sub_into_this(Big_Int const& left, Digits const& right, int right_sign)
        {
            Location loc = (left.slice().data == _digits.data()) 
                ? Location::IN_PLACE 
                : Location::OUT_OF_PLACE;

            set_sign(left.sign());
            if(left.sign() == right_sign)
            {
                grow_buffer(required_add_out_size(left.size(), right.size()));
                Slice my_buf = buffer();
                Slice sum = add<Digit>(&my_buf, left.slice(), big_int::slice(right), loc);
         
                resize_slice(sum.size);
                return;
            }
                
            Batch_Overflow diff;
            if(left.size() >= right.size())
            {
                grow_buffer(required_sub_out_size(left.size(), right.size()));
                Slice my_buf = buffer();
                diff = batch_sub_overflow_long<Digit>(&my_buf, left.slice(), big_int::slice(right), loc);
            }
            else
            {
                //flip args and flip sign
                flip_sign();
                grow_buffer(required_sub_out_size(right.size(), left.size()));
                Slice my_buf = buffer();
                diff = batch_sub_overflow_long<Digit>(&my_buf, big_int::slice(right), left.slice(), Location::OUT_OF_PLACE);
            }

            //if result is negative in total, place the overflown ammount into the digits
            // and flip sign
            if(diff.overflow != 0)
            {
                flip_sign();
                resize_slice(1);
                slice()[0] = 0 - slice()[0];
            }
            else
            {
                Slice stripped = striped_trailing_zeros(diff.slice);
                resize_slice(stripped.size);
            }
        }

    public:
        void add_into_this(Big_Int const& left, Big_Int const& right)
        {
            return add_sub_into_this(left, right._digits, right._sign);
        }

        void sub_into_this(Big_Int const& left, Big_Int const& right)
        {
            return add_sub_into_this(left, right._digits, -right._sign);
        }

        void mul_into_this(Big_Int const& left, Big_Int const& right)
        {
            size_t required_out = required_mul_out_size(left.size(), right.size());
            size_t required_aux = required_mul_karatsuba_aux_size(left.size(), right.size(), optims.max_recursion_depth);
            grow_buffer(required_out);
            grow_aux(required_aux);

            Slice my_buf = buffer();
            Slice my_aux = aux();
            Slice out = mul(&my_buf, &my_aux, left.slice(), right.slice(), optims);

            set_sign(left.sign() * right.sign());
            resize_slice(out.size);
        };

        void div_rem_into_this(Big_Int const& left, Big_Int const& right)
        {
            if(right.size() == 0)
                throw std::domain_error("dvision by 0!");

            size_t required_out = required_div_quotient_size(left.size(), right.size());
            size_t required_aux = required_div_remainder_size(left.size(), right.size());
            grow_buffer(required_out);
            grow_aux(required_aux);

            Slice quotient = buffer();
            Slice remainder = aux();

            Div_Result out = div(&quotient, &remainder, left.slice(), right.slice(), optims);
            set_sign(left.sign() * right.sign());
            resize_slice(out.quotient.size);
            _aux.resize_for_overwrite(out.remainder.size);
        };

        void rem_into_this(Big_Int const& left, Big_Int const& right)
        {
            if(right.size() == 0)
                throw std::domain_error("remainder of 0!");

            size_t required_remainder = required_div_remainder_size(left.size(), right.size());
            grow_buffer(required_remainder);

            Slice my_buff = buffer();
            Slice out = rem(&my_buff, left.slice(), right.slice(), optims);
            set_sign(left.sign());
            resize_slice(out.size);
        };

        void pow_into_this(Big_Int const& num, imax power)
        {
            if(num.sign() == -1 && power % 2 == 0)
                set_sign(1);

            if((size_t) power < 0)
                throw std::domain_error("negative powers not supported!");

            size_t log2_digits = log2(num.slice());
            size_t required_out = required_pow_out_size(log2_digits, (size_t) power, BIT_SIZE<umax>);
            size_t required_aux = optimal_pow_aux_size(log2_digits, (size_t) power, BIT_SIZE<umax>);
            grow_buffer(required_out);
            grow_aux(required_aux);

            Slice my_buf = buffer();
            Slice my_aux = aux();

            Slice res = pow_by_squaring<umax>(&my_buf, &my_aux, num.slice(), (size_t) power, optims);
            resize_slice(res.size);
        }

        void root_into_this(Big_Int const& num, imax power)
        {
            if(num.sign() == -1 && power % 2 == 0)
                throw std::domain_error("even root of negative number!");

            if(power < 0)
                throw std::domain_error("negative root of number!");

            if(tiny_num::high_bits((umax) power) != 0)
                throw std::domain_error("root too big! (numbers bigger than 32 bits currently not supported)");

            size_t log2_digits = log2(num.slice());
            size_t required_out = required_root_out_size(log2_digits, (size_t) power, BIT_SIZE<umax>);
            size_t required_aux = optimal_root_aux_size(log2_digits, (size_t) power, BIT_SIZE<umax>);
            grow_buffer(required_out);
            grow_aux(required_aux);

            Slice my_buf = buffer();
            Slice my_aux = aux();

            Slice res = root<umax>(&my_buf, &my_aux, num.slice(), (size_t) power, optims);
            resize_slice(res.size);
        }

        void fact_into_this(imax of_number)
        {
            if(of_number < 0)
                throw std::domain_error("factorial of negative number!");

            size_t required_out = required_factorial_size((size_t) of_number);
            grow_buffer(required_out);
            Slice my_buf = buffer();
            Slice res = factorial<umax>(&my_buf, (size_t) of_number, optims);
            resize_slice(res.size);
        }

        //get the number as single digit or throw error
        imax to_single_digit()
        {
            if(digit_count() != 1 || (imax) slice()[0] < 0)
                throw std::domain_error("is not convertile to single digit!");

            return slice()[0] * sign();
        }

    public:
        bool operator ==(Big_Int const& other) const
        {
            if(this->_sign != other._sign)
                return false;

            return compare(slice(), other.slice()) == 0;
        }

        bool operator !=(Big_Int const& other) const
        {
            return (*this == other) == false;
        }

        auto operator <=>(Big_Int const& other) const
        {
            //compare(this, other) --> -1 <=> this < other
            //                          0 <=> this = other
            //                          1 <=> this > other

            // -1 ?  1 --> -1
            //  1 ? -1 -->  1
            if(this->_sign != other._sign)
                return this->_sign;

            int res = compare(slice(), other.slice());
            return res * this->_sign;
        }

        Big_Int& operator +=(Big_Int const& other)
        {
            add_into_this(*this, other);
            return *this;
        }

        Big_Int& operator -=(Big_Int const& other)
        {
            sub_into_this(*this, other);
            return *this;
        }

        Big_Int& operator *=(Big_Int const& other)
        {
            mul_into_this(*this, other);
            return *this;
        }

        Big_Int& operator /=(Big_Int const& other)
        {
            div_rem_into_this(*this, other);
            return *this;
        }

        Big_Int& operator %=(Big_Int const& other)
        {
            rem_into_this(*this, other);
            return *this;
        }

        Big_Int operator +(Big_Int const& other) const
        {
            Big_Int out;
            out.add_into_this(*this, other);
            return out;
        }

        Big_Int operator -(Big_Int const& other) const
        {
            Big_Int out;
            out.sub_into_this(*this, other);
            return out;
        }

        Big_Int operator *(Big_Int const& other) const
        {
            Big_Int out;
            out.mul_into_this(*this, other);
            return out;
        }

        Big_Int operator /(Big_Int const& other) const
        {
            Big_Int out;
            out.div_rem_into_this(*this, other);
            return out;
        }

        Big_Int operator %(Big_Int const& other) const
        {
            Big_Int out;
            out.rem_into_this(*this, other);
            return out;
        }

        Big_Int operator -() const
        {
            Big_Int copy = *this;
            copy.flip_sign();
            return copy;
        }

        friend class Big_Int_Popper;
    };

    class Big_Int_Popper
    {
        Pop_Digit<umax> popper;
        Big_Int num;

        public:
        Big_Int_Popper(Big_Int big, size_t base = 10)
        {
            num = std::move(big);
            num.grow_aux(num.size());

            Slice my_buffer = num.buffer();
            Slice my_temp = num.aux();
            popper = pop_digit_init<umax>(&my_buffer, my_temp, base);
        }

        int get_sign() {
            return num.sign();
        }

        bool ok() {
            return popper.ok;
        }

        auto error_type() {
            return popper.error;
        }

        umax pop(){
            return pop_digit(&popper, optims);
        }
    };
}

using big_int::Digit;
using big_int::Big_Int;
using big_int::Big_Int_Pusher;
using big_int::Big_Int_Popper;

Big_Int pow(Big_Int const& left, big_int::imax power)
{
    Big_Int out;
    out.pow_into_this(left, power);
    return out;
}

Big_Int root(Big_Int const& left, big_int::imax power)
{
    Big_Int out;
    out.root_into_this(left, power);
    return out;
}

Big_Int sqrt(Big_Int const& left)
{
    return root(left, 2);
}

Big_Int factorial(big_int::imax of_number)
{
    Big_Int out;
    out.fact_into_this(of_number);
    return out;
}
namespace big_int
{
    void test_big_int()
    {
        std::cout << "TESTING BIG INT "<< std::endl;
        Big_Int a = factorial(17);
        Big_Int b = factorial(25);

        assert(a - a == 0);
        assert(a + (-a) == 0);
        assert((-a) - a == a*-2);
        assert(a - b != a);
        assert(a - b == -(b - a));
        assert(-a - b == -(b + a));

        assert(a*2 == a + a);
        assert(a*2 == a + a);
        assert(a*2 == a + a);
        assert((a + b)*(a - b) == pow(a, 2) - pow(b, 2));
        assert(pow((a + b), 2) == (a*a + a*b*2 + b*b));

        assert(root(pow(b, 3), 3) == b);
        std::cout << "=== OK ===\n" << std::endl;      
    }
}