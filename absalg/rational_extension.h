#pragma once
#include "field.h"
#include <iostream>
#include <real_field.h>

template <typename R>
class rational_extension: public field
{
public:
	rational_extension(R _p) :rational_extension(_p, R::_1()) {}
	rational_extension(R _p, R _q) :p(_p), q(_q) { reduce(); }
	rational_extension(int _p=0,int _q=1):p(_p),q(_q){}

	bool operator!=(const rational_extension& s) const
	{
		return (s.p != p) || (s.q != q);
	}
	bool operator==(const rational_extension& s) const
	{
		return !(*this != s);
	}
	inline static rational_extension _0()
	{
		return rational_extension();
	}
	inline static rational_extension _1()
	{
		return rational_extension(R::_1());
	}
	const R& nominator() const
	{
		return p;
	}

	const R& denominator() const
	{
		return q;
	}
	bool is_zero() const
	{
		return p.is_zero();
	}

	bool is_one() const
	{
		return p == q;
	}
	explicit operator R()
	{
		return p.div(q);
	}
	rational_extension& operator+=(const rational_extension& a)
	{
		p = p * a.q + a.p * q;
		q *= a.q;
		reduce();
		return *this;
	}
	rational_extension& operator-=(const rational_extension& a)
	{
		p = p * a.q - a.p * q;
		q *= a.q;
		reduce();
		return *this;
	}

	rational_extension operator-() const
	{
		return rational_extension(-p, q);
	}

	rational_extension& operator*=(const rational_extension& a)
	{
		p *= a.p;
		q *=a.q;
		reduce();
		return *this;
	}
	rational_extension& operator/=(const rational_extension& a)
	{
		p *= a.q;
		q *=a.p;
		reduce();
		return *this;
	}
	rational_extension& operator+=(const R& a)
	{
		return this->operator+=(rational_extension(a));
	}
	rational_extension& operator-=(const R& a)
	{
		return this->operator-=(rational_extension(a));
	}

	rational_extension& operator*=(const R & a)
	{
		p *= a;
		reduce();
		return *this;
	}
	rational_extension& operator/=(const R & a)
	{
		q*=a;
		reduce();
		return *this;
	}

	rational_extension& operator+=(int a)
	{
		return this->operator+=(R(a));
	}
	rational_extension& operator-=(int a)
	{
		return this->operator-=(R(a));
	}

	rational_extension& operator*=(int a)
	{
		p *= a;
		reduce();
		return *this;
	}
	rational_extension& operator/=(int a)
	{
		q *= a;
		reduce();
		return *this;
	}


	real_field abs() const requires ring_constraints::has_abs<R> 
	{
		return p.abs() / q.abs();
	}

private:
	void reduce()
	{
		if (p == R::_0())
		{
			q = R::_1();
			return;
		}
		R d(R::gcd(p, q));
		p = p.div(d);
		q = q.div(d);
	}
	R p, q;
};

template<typename R>
rational_extension<R> operator+(const rational_extension<R>& a, const rational_extension<R>& b)
{
	rational_extension c(a);
	return c += b;
}

template<typename R>
rational_extension<R> operator-(const rational_extension<R>& a, const rational_extension<R>& b)
{
	rational_extension c(a);
	return c -= b;
}

template<typename R>
rational_extension<R> operator*(const rational_extension<R>& a, const rational_extension<R>& b)
{
	rational_extension c(a);
	return c *= b;
}

template<typename R>
rational_extension<R> operator/(const rational_extension<R>& a, const rational_extension<R>& b)
{
	rational_extension c(a);
	return c /= b;
}

template<typename R>
std::ostream& operator<<(std::ostream& H, const rational_extension<R> a)
{
	return H << "[" << a.nominator() << "|" << a.denominator() << "]";
}