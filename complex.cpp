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

complex math_rz::complex::inner_product(const complex& z) const
{
    return (*this) * z.conj();
}

real_field math_rz::complex::norm() const
{
    return abs();
}

complex& math_rz::complex::operator/=(const real_field& s)
{
    real() /= s;
    imag() /= s;
    return *this;
}

/*complex math_rz::operator/(const complex &s,const real_field& a)
{
    complex R(s);
    return R /= a;
}*/

complex math_rz::complex::operator-() const
{
    return complex(-real(), -imag());
}

complex& math_rz::complex::operator*=(int a)
{
    this->_Val[0] *= a;
    this->_Val[1] *= a;
    return  *this;
}


complex math_rz::operator""_c(long double a)
{
    return complex(a,0);
}

complex math_rz::operator""_c(unsigned long long a)
{
    return complex(a,0);
}

complex math_rz::operator*(int a, const complex& s)
{
    auto w = s;
    return w *= a;
}

complex math_rz::operator*(const complex& s,int a)
{
    return a * s;
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

complex& complex::operator*=(const std::complex<real_field>& s)
{
    return *this = (*this) * s;
}

complex& complex::operator/=(const std::complex<real_field>& s)
{
    return *this = (*this) / s;
}