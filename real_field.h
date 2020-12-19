#pragma once
#pragma once
#include "absalg/field.h"
#include <algorithm>
#include <iostream>
#include <compare>

namespace math_rz
{
    class real_field :
        public field
    {
    public:
        real_field(long double m);
        real_field(double m);
        real_field(float m);
        real_field(int m = 0);
        real_field(unsigned long long m);
        inline static const real_field _0() { return real_field(); };
        inline static const real_field _1() { return real_field(1); };
        real_field mod(const integral_ring& a) const {
            return 0;
        }

        real_field abs() const;
        bool is_zero() const;
        bool is_one() const;

        real_field div(const integral_ring& s) { return 0; }
        const real_field& I0() const
        {
            return _0();
        }
        const real_field& I1() const
        {
            return _1();
        }


        operator long double& () { return _v; };
        operator long const double& ()const { return _v; };
        long double _v;
    };

}