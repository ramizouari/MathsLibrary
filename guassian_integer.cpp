#include "guassian_integer.h"

guassian_integer::guassian_integer(const integer& a, const integer& b) :std::complex<integer>(a, b)
{
}

bool guassian_integer::is_zero() const
{
    return real().is_zero()&&imag().is_zero();
}

bool guassian_integer::is_one() const
{
    return real().is_one()&&imag().is_zero();
}
