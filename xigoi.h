#pragma once

#include <algorithm>
#include <cassert>
#include <exception>
#include <execution>
#include <iostream>
#include <ranges>
#include <vector>

namespace xigoi
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
    using naf::Iota_View;
    using std::views::reverse;

    using Digit = long long;
    using Digits = std::vector<Digit>;
    using Size = Digits::size_type;

    const auto parallel = std::execution::par;

    /// An arbitrary-size integer, stored using a redundant balanced radix system.
    class Integer {

    private:
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
            Size length = max(digits.size(), other.digits.size());
            digits.resize(length + 1);
            transform(parallel, other.digits.begin(), other.digits.end(),
                digits.begin(), digits.begin(),
                [](Digit x, Digit y) { return x + y; });
            Digits carry(length);
            
            
            auto indices = Iota_View(Size(0), length);
            for_each(parallel, indices.begin(), indices.end(), [&](Size index) {
                if (digits[index] >= max_digit) {
                    carry[index] = 1;
                    digits[index] -= base;
                } else if (digits[index] <= -max_digit) {
                    carry[index] = -1;
                    digits[index] += base;
                } else {
                    carry[index] = 0;
                }
            });

            transform(parallel, carry.begin(), carry.end(), digits.begin() + 1,
                digits.begin() + 1, [](Digit x, Digit y) { return x + y; });

                while (digits.empty() == false && digits.back() == 0) {
                    digits.pop_back();
                }

            return *this;
        }

        /// Add two Integers.
        Integer operator+(const Integer &other) const {
            auto result = *this;
            result += other;
            return result;
        }

        /// Calculate the additive inverse of an Integer in place.
        void negate() {
            transform(parallel, digits.begin(), digits.end(), digits.begin(),
                [](Digit digit) { return -digit; });
        }

        /// Calculate the additive inverse of an Integer.
        Integer operator-() const {
            Digits result(digits.size());
            transform(parallel, digits.begin(), digits.end(), result.begin(),
                [](Digit digit) { return -digit; });
            return Integer(result);
        }

        /// Subtract two Integers in place.
        Integer &operator-=(const Integer &other) {
            *this += (-other);
            return *this;
        }

        /// Subtract two Integers.
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

        /// Multiply two Integers.
        Integer operator*(const Integer &other) const {
            // Optimization: this algorithm is faster
            // when multiplying a shorter number by a longer number.
            if (digits.size() > other.digits.size()) {
                return other * *this;
            }
            Size result_length = digits.size() + other.digits.size();
            Integer result = 0;
            Digits augend(result_length);
            Digits carry(result_length);
            auto this_indices = Iota_View(Size(0), digits.size());
            auto other_indices = Iota_View(Size(0), other.digits.size());
            for_each(this_indices.begin(), this_indices.end(), [&](Size i) {
                fill(parallel, augend.begin(), augend.end(), 0);
                fill(parallel, carry.begin(), carry.end(), 0);
                for_each(parallel, other_indices.begin(), other_indices.end(),
                    [&](Size j) {
                        auto product = digits[i] * other.digits[j];
                        augend[i + j] = product % base;
                        carry[i + j + 1] = product / base;
                    });
                result += Integer(augend);
                result += Integer(carry);
                });
            return result;
        }

        /// Multiply two Integers in place.
        Integer &operator*=(const Integer &other) {
            *this = *this * other;
            return *this;
        }

        /// Divide an Integer by two, rounding toward zero.
        Integer half() const {
            //@BUG - when zero input acesses ouyt of bounds
            if(this->digits.size() == 0)
                return Integer();

            Integer result(Digits(digits.size()));
            transform(parallel, digits.begin(), digits.end(), result.digits.begin(),
                [](Digit digit) { return digit / 2; });
            transform(parallel, digits.begin() + 1, digits.end(), result.digits.begin(),
                result.digits.begin(), [](Digit high, Digit low) {
                    return low + (high % 2) * half_base;
                });
            if (result.digits.back() == 0) {
                result.digits.pop_back();
            }
            return result;
        }

        /// Divide two Integers, rounding toward zero.
        /// Uses binary search because normal division
        /// would be extremely complicated and
        /// this works quite well.
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
            Integer old_low;
            Integer old_high;
            Integer diff;

            do 
            {
                middle = (low + high).half();
                diff = other * middle - *this;

                if (diff.is_zero()) {
                    return middle;
                } else if (diff.is_positive()) {
                    high = middle;
                } else {
                    low = middle + 1;
                }
            } 
            while (high != low);
            return low - 1;
        }

        /// Divide two Integers in place, rounding toward zero.
        Integer &operator/=(const Integer &other) {
            *this = *this / other;
            return *this;
        }

        /// Calculate the remainder after division of two Integers.
        /// m / n + m % n == m
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
        auto range = Iota_View(Digit(1), n + 1);
        return accumulate(range.begin(), range.end(), Integer(1),
            [](const Integer &i, Digit k) { return i * k; });
    }

    // As a small test, calculate and print the factorial of 1000.
    // Takes about 10 seconds on my machine, with most of the time
    // spent converting to decimal.
    int main() { cout << factorial(1000) << "\n"; return 0; }
}