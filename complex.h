#pragma once
#include <complex>
#include "field.h"
#include "real_field.h"
class complex:public std::complex<real_field>,public field
{
public:
	complex(const std::complex<real_field>& a) :std::complex<real_field>(a) {}
	complex(const real_field& a = 0, const real_field& b = 0);
	complex(int a);
	static complex _0() { return std::complex<real_field >(); }
	static complex _1() { return std::complex<real_field>(1); }
	bool is_zero() const;
	bool is_one() const;
	inline operator std::complex<real_field>() { return *this; }
};

::complex operator""_c(long double a);
::complex operator""_c(unsigned long long a);
