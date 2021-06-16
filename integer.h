#pragma once
#include "absalg/integral_ring.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <compare>
namespace math_rz {
    class real_field;
    class integer 
    {
    public:
        constexpr integer(const int m = 0):_v(m) {};
        constexpr integer(long long int m) :_v(m) {};
        constexpr integer(const double s) :_v(s) {};
        inline static integer _0() { return integer(); };
        inline static integer _1() { return integer(1); };
        integer mod(const integer& a) const;
        integer div(const integer& a) const;
        inline static integer gcd(const integer& a, const integer& b) { return std::gcd(a._v, b._v); }
        inline static integer lcm(const integer& a, const integer& b) { return std::lcm(a._v, b._v); }
        inline static std::pair<integer, integer> euclidean_division(const integer& a, const integer& b) {
            auto R(std::div(a._v, b._v));
            return std::make_pair<integer, integer>(R.quot, R.rem);
        }

        real_field abs() const;
        real_field inner_product(const real_field& other) const;
        inline integer operator-() const { return -_v; }
        integer& operator+=(const integer& o);

        inline operator long long& () { return _v; }
        inline operator const long long& () const { return _v; }
        inline explicit operator double() { return _v; }
        inline explicit operator long double() { return _v; }
        inline explicit operator float() { return _v; }
        bool is_zero() const;
        bool is_one() const;
        long long _v;
    };
    //std::ostream& operator<<(std::ostream& H, const integer& a);
}