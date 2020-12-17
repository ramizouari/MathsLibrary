#include "real_field.h"

real_field::real_field(double m) :_v(m)
{
}

bool real_field::is_zero() const
{
	return _v==0;
}

bool real_field::is_one() const
{
	return _v==1;
}

bool real_field::operator==(const real_field& a) const
{
	return _v == a._v;
}

bool real_field::operator!=(const real_field& a) const
{
	return _v != a._v;
}

real_field real_field::operator-()
{
	return -_v;
}

real_field& real_field::operator+()
{
	return *this;
}

real_field& real_field::operator+=(const real_field& a)
{
	_v += a._v;
	return *this;
}

real_field& real_field::operator-=(const real_field& a)
{
	_v -= a._v;
	return *this;
}

real_field& real_field::operator*=(const real_field& a)
{
	_v *= a._v;
	return *this;
}

real_field& real_field::operator/=(const real_field& a)
{
	_v /= a._v;
	return *this;
}

real_field& real_field::operator%=(const real_field& a)
{
	_v = 0;
	return *this;
}

std::partial_ordering real_field::operator<=>(const real_field& a)
{
	return _v <=> a._v;
}

ring& real_field::operator+=(int n)
{
	_v += n;
	return *this;
}

ring& real_field::operator-=(int n)
{
	_v -= n;
	return *this;
}

ring& real_field::operator*=(int n)
{
	_v *= n;
	return *this;
}

real_field operator+(const real_field& a, const real_field& b)
{
	real_field c = a;
	return c += b;
}

real_field operator-(const real_field& a, const real_field& b)
{
	return a._v - b._v;
}

real_field operator*(const real_field& a, const real_field& b)
{
	real_field c = a;
	return c *= b;
}

real_field operator/(const real_field& a, const real_field& b)
{
	real_field c = a;
	return c /= b;
}

real_field operator%(const real_field& a, const real_field& b)
{
	real_field c = a;
	return c %= b;
}

std::ostream& operator<<(std::ostream& H, const real_field& a)
{
	return H << a._v;
}