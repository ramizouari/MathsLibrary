#pragma once
#include "integral_ring.h"
#include <algorithm>
#include <iostream>
class integer :
    public integral_ring
{
public:
    integer(const int m = 0);
    inline static const integer _0() { return integer(); };
    inline static const integer _1() { return integer(1); };
    integer&& mod(const integral_ring& a) {
        const integer& b = dynamic_cast<const integer&>(a);
        return integer(this->_v % b._v);
    }
    integer&& div(const integral_ring&) { return integer(0); }
    const integer& I0() const
    {
        return _0();
    }
    const integer& I1() const
    {
        return _1();
    }

    bool operator==(const integer& a) const;
    operator int&() { return _v; };
    integer operator-();
    integer& operator+=(const integer& a);
    integer& operator-=(const integer& a);
    integer& operator*=(const integer& a);
    integer& operator/=(const integer& a);
    integer& operator%=(const integer& a);
    virtual ring& operator+=(int n);
    virtual ring& operator-=(int n);
    virtual ring& operator*=(int n);
    int _v;
};

integer operator+(const integer& a, const integer& b);
integer operator*(const integer& a, const integer& b);
integer operator/(const integer& a, const integer& b);
integer operator%(const integer& a, const integer& b);
std::ostream& operator<<(std::ostream& H, const integer& a);