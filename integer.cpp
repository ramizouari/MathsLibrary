#include "integer.h"

integer::integer(const int m):_v(m)
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

bool integer::operator==(const integer& a) const
{
	return _v==a._v;
}

bool integer::is_zero() const
{
	return _v==0;
}

bool integer::is_one() const
{
	return _v == 1;;
}

integer integer::operator-() const
{
	return -_v;
}

const integer& integer::operator+() const
{
	return *this;
}

integer& integer::operator+=(const integer& a)
{
	_v += a._v;
	return *this;
}

integer& integer::operator-=(const integer& a)
{
	_v -= a._v;
	return *this;
}

integer& integer::operator*=(const integer& a)
{
	_v *= a._v;
	return *this;
}

integer& integer::operator/=(const integer& a)
{
	_v /= a._v;
	return *this;
}

integer& integer::operator%=(const integer& a)
{
	_v %= a._v;
	return *this;
}

std::strong_ordering integer::operator<=>(const integer& b)
{
	return _v <=> b._v;
}

ring& integer::operator+=(int n)
{
	_v += n;
	return *this;
}

ring& integer::operator-=(int n)
{
	_v -= n;
	return *this;
}

ring& integer::operator*=(int n)
{
	_v *= n;
	return *this;
}

integer operator+(const integer& a, const integer& b)
{
	integer c = a;
	return c += b;
}

integer operator-(const integer& a, const integer& b)
{
	return a._v - b._v;
}

integer operator*(const integer& a, const integer& b)
{
	integer c = a;
	return c *= b;
}

integer operator/(const integer& a, const integer& b)
{
	integer c = a;
	return c /= b;
}

integer operator%(const integer& a, const integer& b)
{
	integer c = a;
	return c %= b;
}

std::ostream& operator<<(std::ostream& H, const integer& a)
{
	return H << a._v;
}