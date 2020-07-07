#include "integer.h"

integer::integer(const int m):_v(m)
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

integer integer::operator-()
{
	return -_v;
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