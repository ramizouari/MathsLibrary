#include "complex.h"


using namespace math_rz;
complex::complex(const real_field& a, const real_field& b):std::complex<real_field>(a,b)
{
}

math_rz::complex::complex(double s):std::complex<real_field>(s)
{
}

math_rz::complex::complex(float s): std::complex<real_field>(s)
{
}

math_rz::complex::complex(long double s): std::complex<real_field>(s)
{
}

complex::complex(int a) : std::complex<real_field>(a)
{
}

bool complex::is_zero() const
{
    return real().is_zero() && imag().is_zero();
}

bool complex::is_one() const
{
    return real().is_one() && imag().is_zero();
}

real_field complex::abs() const
{
    return std::abs(*this);
}

complex complex::conj() const
{
    return std::conj(*this);
}

complex operator""_c(long double a)
{
    return complex(a);
}

complex operator""_c(unsigned long long a)
{
    return complex(a,0);
}

complex operator""_i(long double a)
{
    return complex(0,a);
}
