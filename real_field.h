#pragma once
#include "absalg/field.h"
#include <algorithm>
#include <iostream>
#include <compare>
#include "integer.h"

namespace math_rz
{
    class real_field :
        public field
    {
    public:
        inline static constexpr int dimension = 1;
        using base_field = real_field;
        real_field(long double m);
        real_field(double m);
        real_field(float m);
        real_field(int m = 0);
        real_field(unsigned long long m);
        real_field(long long m);
        real_field(const integer& a);
        inline static const real_field _0() { return real_field(); };
        inline static const real_field _1() { return real_field(1); };
        static bool exact;
        real_field mod(const real_field& a) const {
            return 0;
        }

        real_field abs() const;
        real_field norm() const;
        real_field conj() const;
        real_field inner_product(const real_field& a);
        real_field inv()const;
        bool is_zero() const;
        bool is_one() const;

        real_field div(const real_field& s) { return _v/ s._v; }
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