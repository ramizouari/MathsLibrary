#include "guassian_integer.h"

using namespace math_rz;
guassian_integer::guassian_integer(const integer& a, const integer& b) :std::complex<integer>(a, b)
{
}

guassian_integer::guassian_integer(int a):std::complex<integer>(a)
{
}

math_rz::guassian_integer::guassian_integer(const guassian_integer& a, int b):guassian_integer(b)
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

guassian_integer guassian_integer::conj() const
{
    return guassian_integer(real(),-imag());
}
