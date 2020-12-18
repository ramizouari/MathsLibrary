#include "complex.h"

complex::complex(const real_field& a, const real_field& b):std::complex<real_field>(a,b)
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

::complex operator""_c(long double a)
{
    return ::complex(a);
}

::complex operator""_c(unsigned long long a)
{
    return ::complex(a,0);
}

::complex operator""_i(long double a)
{
    return ::complex(0,a);
}
