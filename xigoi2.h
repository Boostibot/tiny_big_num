#pragma once

#include <algorithm>
#include <cassert>
#include <exception>
#include <execution>
#include <iostream>
#include <ranges>
#include <vector>

#include "naf.h"

namespace xigoi2
{
    using std::accumulate;
    using std::cout;
    using std::domain_error;
    using std::fill;
    using std::for_each;
    using std::make_pair;
    using std::max;
    using std::ostream;
    using std::out_of_range;
    using std::pair;
    using std::string;
    using std::transform;
    using std::ranges::iota_view;
    using std::views::reverse;

    using Digit = naf::Digit;
    using Digits = naf::Digits;
    using Size = naf::Size;

    const auto parallel = std::execution::par;

    /// An arbitrary-size integer, stored using a redundant balanced radix system.
    class Integer {

    public:
        /// Digits are stored starting with the least significant position.
        Digits digits;
        /// The radix used for the integers.
        static const Digit base = Digit(1) << 31;
        static const Digit half_base = Digit(1) << 30;
        /// The individual digits can range from -max_digit to max_digit inclusive.
        static const Digit max_digit = base - 1;

        /// If the Integer is equal to a single-digit Integer,
        /// returns the single digit.
        /// If not, the behavior is undefined!
        Digit unsafe_to_digit() const {
            switch (digits.size()) {
            case 0:
                return 0;
            case 1:
                return digits[0];
            default:
                if (digits[0] < 0) {
                    return digits[0] + base;
                } else {
                    return digits[0] - base;
                }
            }
        }

    public:
        /// Empty constructor
        Integer() {}

        /// Convert a single Digit to an Integer.
        Integer(Digit x) { // cppcheck-suppress noExplicitConstructor
            while (x != 0) {
                digits.push_back(x % base);
                x /= base;
            }
        }

        /// Convert a list of Digits to an Integer.
        explicit Integer(const Digits &digits_) : digits(digits_) {}

        /// Add two Integers in place.
        /// Uses Avizienis' parallel addition algorithm, as described on page 4 in
        /// https://arxiv.org/pdf/1102.5683.pdf
        Integer &operator+=(const Integer &other) {
            Digits carry;
            naf::incr_or_decr(&carry, &this->digits, other.digits, true);

            return *this;
        }

        /// Add two Integers.
        Integer operator+(const Integer &other) const {
            Integer out;
            Digits carry;
            naf::add_or_sub(&out.digits, &carry, this->digits, other.digits, true);
            return out;
            //auto result = *this;
            //result += other;
            //return result;
        }

        /// Calculate the additive inverse of an Integer in place.
        void negate() {
            transform(parallel, digits.begin(), digits.end(), digits.begin(),
                [](Digit digit) { return -digit; });
        }

        /// Calculate the additive inverse of an Integer.
        Integer operator-() const {
            Digits result;
            
            result.resize(digits.size());
            transform(parallel, digits.begin(), digits.end(), result.begin(),
                [](Digit digit) { return -digit; });
            return Integer(result);
        }

        /// Subtract two Integers in place.
        Integer &operator-=(const Integer &other) {
            Digits carry;
            naf::incr_or_decr(&carry, &this->digits, other.digits, false);

            return *this;

            //*this += (-other);
            //return *this;
        }

        /*
        Integer operator-(const Integer &other) const {
            Integer out;
            Digits carry;
            naf::add_or_sub(&out.digits, &carry, this->digits, other.digits, false);
            return out;
        }
        */

         //Subtract two Integers.
        Integer operator-(const Integer &other) const { return *this + (-other); }

        /// Compare an Integer with zero.
        bool is_zero() const { return digits.size() == 0; }
        bool is_nonzero() const { return digits.size() != 0; }
        bool is_positive() const { return digits.size() != 0 and digits.back() > 0; }
        bool is_negative() const { return digits.size() != 0 and digits.back() < 0; }
        bool is_zero_or_positive() const {
            return digits.size() == 0 or digits.back() > 0;
        }
        bool is_zero_or_negative() const {
            return digits.size() == 0 or digits.back() < 0;
        }

        /// Compare two Integers.
        bool operator==(const Integer &other) const {
            return (*this - other).is_zero();
        }
        bool operator!=(const Integer &other) const {
            return (*this - other).is_nonzero();
        }
        bool operator>(const Integer &other) const {
            return (*this - other).is_positive();
        }
        bool operator<(const Integer &other) const {
            return (*this - other).is_negative();
        }
        bool operator>=(const Integer &other) const {
            return (*this - other).is_zero_or_positive();
        }
        bool operator<=(const Integer &other) const {
            return (*this - other).is_zero_or_negative();
        }
        
        Integer operator*(const Integer &other) const {
            Integer out;
            Digits temp;
            Digits carry;
            Digits augend;
            naf::mul(&out.digits, &temp, &augend, &carry, this->digits, other.digits);

            return out;
        }

        /// Multiply two Integers in place.
        Integer &operator*=(const Integer &other) {
            *this = *this * other;
            return *this;
        }

        Integer half() const 
        {
            Integer result;
            naf::half(&result.digits, this->digits);
            return result;
        }

        /// Divide two Integers, rounding toward zero.
        /// Uses binary search because normal division
        /// would be extremely complicated and
        /// this works quite well.
        
        /*
        Integer operator/(const Integer &other) const {
            if (other.is_zero()) {
                throw domain_error("Division by zero");
            }
            if (other.is_negative()) {
                if (is_negative()) {
                    return (-*this) / (-other);
                } else {
                    return -(*this / (-other));
                }
            } else if (is_negative()) {
                return -((-*this) / other);
            }
            Integer low = 0;
            Integer high = *this;
            Integer middle;
            Integer diff;
            do {
                middle = (low + high).half();
                diff = other * middle - *this;
                if (diff.is_zero()) {
                    return middle;
                } else if (diff.is_positive()) {
                    high = middle;
                } else {
                    low = middle + 1;
                }
            } while (low != high);
            return low - 1;
        }
        */

        pair<Integer, Integer> divmod(const Integer &other, bool only_div = false) const
        {
            if (other.is_zero()) {
                throw domain_error("Division by zero");
            }

            Digit sign = 1;
            if (other.is_negative()) {
                if (is_negative()) {
                    sign = 1;
                } else {
                    sign = -1;
                }
            } else if (is_negative()) {
                sign = -1;
            }

            Integer normalized_this = this->is_zero_or_positive() ? *this : -*this;
            Integer normalized_other = other.is_zero_or_positive() ? other : -other;

            Integer diff;
            Integer low; 
            Integer high = normalized_this;
            Integer old_low;
            Integer old_high;

            Integer middle; 
            Integer muled;
            Integer temp;
            Integer augend; 
            Integer carry; 
            Integer one = 1;

            //Integer low = 0;
            //Integer high = *this;
            //Integer middle;
            //Integer diff;
            while(true) 
            {
                //middle = (low + high).half();
                naf::add_or_sub(&temp.digits, &carry.digits, low.digits, high.digits, true);
                naf::half(&middle.digits, temp.digits);

                //diff = other * middle - *this;
                naf::mul(&muled.digits, &temp.digits, &augend.digits, &carry.digits, normalized_other.digits, middle.digits);
                naf::add_or_sub(&diff.digits, &carry.digits, muled.digits, normalized_this.digits, false);

                if (diff.is_zero()) 
                {
                   return make_pair(std::move(middle), std::move(diff));

                } else if (diff.is_positive()) {
                    high = middle;
                } else {
                    naf::add_or_sub(&low.digits, &carry.digits, middle.digits, one.digits, true);
                }
            } 
            while (naf::are_equal(low.digits, high.digits) == false);


            naf::incr_or_decr(&carry.digits, &low.digits, one.digits, false);
            if(sign == -1)
                naf::negate(&low.digits);

            if(only_div)
                return make_pair(std::move(low), Integer());

            naf::mul(&muled.digits, &temp.digits, &augend.digits, &carry.digits, other.digits, low.digits);
            naf::add_or_sub(&diff.digits, &carry.digits, this->digits, muled.digits, false);

            return make_pair(std::move(low), std::move(diff));
        }


        #if 0
        /// Divide two Integers in place, rounding toward zero.
        Integer &operator/=(const Integer &other) {
            *this = *this / other;
            return *this;
        }

        /// Calculate the remainder after division of two Integers.
        /// m / n + m % n == m
        // -m / n + (-m) % n = -m

        Integer operator%(const Integer &other) {
            Integer quot = *this / other;
            return *this - other * quot;
        }

        /// Calculate the remainder after division of two Integers in place.
        Integer &operator%=(const Integer &other) {
            *this = *this % other;
            return *this;
        }

        /// Divide two integers, rounding toward zero.
        /// Returns both the quotient and the remainder.
        pair<Integer, Integer> divmod(const Integer &other) {
            Integer quot = *this / other;
            return make_pair(quot, *this - other * quot);
        }
        #endif

        /// Output the decimal representation of an Integer to a stream.
        friend ostream &operator<<(ostream &out, Integer n);
    };

    /// Output the decimal representation of an Integer to a stream.
    ostream &operator<<(ostream &out, Integer n) {
        if (n.is_zero()) {
            out << '0';
            return out;
        }
        if (n.is_negative()) {
            out << '-';
            n.negate();
        }
        Digits output_parts;
        while (n.is_nonzero()) {
            auto divmod = n.divmod(1'000'000'000);
            n = divmod.first;
            output_parts.push_back(divmod.second.unsafe_to_digit());
        }
        for (Digit part : reverse(output_parts)) {
            out << part;
            out.width(9);
            out.fill('0');
        }
        out.width(0);
        out.fill(' ');
        return out;
    }

    /// Calculate the factorial of a Digit.
    Integer factorial(Digit n) {
        auto range = iota_view(Digit(1), n + 1);
        return accumulate(range.begin(), range.end(), Integer(1),
            [](const Integer &i, Digit k) { return i * k; });
    }

    // As a small test, calculate and print the factorial of 1000.
    // Takes about 10 seconds on my machine, with most of the time
    // spent converting to decimal.
    int main() { cout << factorial(1000) << "\n"; return 0; }
}