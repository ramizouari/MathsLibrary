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
