#include "integer.h"

integer::integer(const int m):v(m)
{
}

bool integer::operator==(const integer& a) const
{
	return v==a.v;
}

integer integer::operator-()
{
	return -v;
}

integer& integer::operator+=(const integer& a)
{
	v += a.v;
	return *this;
}

integer& integer::operator-=(const integer& a)
{
	v -= a.v;
	return *this;
}

integer& integer::operator*=(const integer& a)
{
	v *= a.v;
	return *this;
}

integer& integer::operator/=(const integer& a)
{
	v /= a.v;
	return *this;
}

integer& integer::operator%=(const integer& a)
{
	v %= a.v;
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
