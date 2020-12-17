#pragma once
#include "integral_ring.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <compare>
class integer :
    public integral_ring
{
public:
    integer(const int m = 0);
    integer(const double s);
    inline static integer _0() { return integer(); };
    inline static integer _1() { return integer(1); };
    integer mod(const integer& a) const;
    integer div(const integer& a) const;
    inline static integer gcd(const integer& a, const integer& b) { return std::gcd(a._v, b._v); }
    inline static integer lcm(const integer& a, const integer& b) { return std::lcm(a._v, b._v); }
    inline static std::pair<integer, integer> euclidean_division(const integer& a, const integer& b) {
        auto R(std::div(a._v, b._v));
        return std::make_pair(R.quot, R.rem);
    }
    bool operator==(const integer& a) const;
    inline explicit operator int&() { return _v; }
    inline explicit operator double() { return _v; }
    inline explicit operator long double() { return _v; }
    inline explicit operator float() { return _v; }
    bool is_zero() const;
    bool is_one() const;
    integer operator-();
    integer& operator+();
    integer& operator+=(const integer& a);
    integer& operator-=(const integer& a);
    integer& operator*=(const integer& a);
    integer& operator/=(const integer& a);
    integer& operator%=(const integer& a);
    std::strong_ordering operator<=>(const integer& b);
    virtual ring& operator+=(int n);
    virtual ring& operator-=(int n);
    virtual ring& operator*=(int n);
    int _v;
};

integer operator+(const integer& a, const integer& b);
integer operator-(const integer& a, const integer& b);
integer operator*(const integer& a, const integer& b);
integer operator/(const integer& a, const integer& b);
integer operator%(const integer& a, const integer& b);
std::ostream& operator<<(std::ostream& H, const integer& a);