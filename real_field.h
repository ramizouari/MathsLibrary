#pragma once
#pragma once
#include "field.h"
#include <algorithm>
#include <iostream>
class real_field :
    public field
{
public:
    real_field(double m = 0);
    inline static const real_field _0() { return real_field(); };
    inline static const real_field _1() { return real_field(1); };
    real_field&& mod(const integral_ring& a) {
        return 0;
    }
    real_field&& div(const integral_ring&) { return real_field(0); }
    const real_field& I0() const
    {
        return _0();
    }
    const real_field& I1() const
    {
        return _1();
    }

    bool operator==(const real_field& a) const;
    bool operator!=(const real_field& a) const;
    operator double& () { return _v; };
    real_field operator-();
    real_field& operator+=(const real_field& a);
    real_field& operator-=(const real_field& a);
    real_field& operator*=(const real_field& a);
    real_field& operator/=(const real_field& a);
    real_field& operator%=(const real_field& a);
    virtual ring& operator+=(int n);
    virtual ring& operator-=(int n);
    virtual ring& operator*=(int n);
    double _v;
};

real_field operator+(const real_field& a, const real_field& b);
real_field operator*(const real_field& a, const real_field& b);
real_field operator/(const real_field& a, const real_field& b);
real_field operator%(const real_field& a, const real_field& b);
std::ostream& operator<<(std::ostream& H, const real_field& a);