#pragma once
#include "absalg/field.h"
#include <algorithm>
#include <iostream>
#include <compare>
#include <functional>
#include "integer.h"

namespace math_rz
{
    class real_field
    {
    public:
        inline static constexpr int dimension = 1;
        using base_field = real_field;
        constexpr real_field(long double m):_v(m) {}
        constexpr real_field(double m) : _v(m) {}
        constexpr real_field(float m) : _v(m) {}
        constexpr real_field(int m = 0) : _v(m) {}
        constexpr real_field(unsigned long long m) : _v(m) {}
        constexpr real_field(long long m) : _v(m) {}
        inline real_field(const integer& a) : _v(a) {}
 
        static bool exact;
        real_field mod(const real_field& a) const {
            return 0;
        }

        real_field abs() const;
        real_field norm() const;
        real_field conj() const;
        real_field inner_product(const real_field& a)const;
        real_field dot_product(const real_field& a)const;
        real_field distance(const real_field &a)const;
        real_field metric(const real_field& a)const;
        real_field inv()const;
        bool is_zero() const;
        bool is_one() const;

        void foreach(const std::function<void(real_field&)>& f);

        real_field div(const real_field& s) { return _v/ s._v; }
        inline constexpr operator long double& () { return _v; };
        inline constexpr operator long const double& ()const { return _v; };
        inline explicit constexpr operator long double ()const { return _v; };
        long double _v;
    };
}