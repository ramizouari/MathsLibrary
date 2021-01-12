#include "integer.h"
#include "real_field.h"


using namespace math_rz;
integer::integer(const int m):_v(m)
{
}

integer::integer(long long int m):_v(m)
{
}

integer::integer(const double s):_v(s)
{
}

integer integer::mod(const integer& a) const
{
	return _v % a._v;
}

integer integer::div(const integer&a) const
{
	return _v/a._v;
}

real_field integer::abs() const
{
	return (long double)std::abs(_v);
}

real_field math_rz::integer::inner_product(const real_field& other) const
{
	return _v * other._v;
}

bool integer::is_zero() const
{
	return _v==0;
}

bool integer::is_one() const
{
	return _v == 1;;
}

std::ostream& operator<<(std::ostream& H, const integer& a)
{
	return H << a._v;
}