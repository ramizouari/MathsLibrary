#pragma once
#include "absalg/ring.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <execution>

template<typename R>
class free_algebra: virtual public ring
{
public:
	free_algebra(){}
	free_algebra(R m):a(1,m){}
	free_algebra(std::vector<R>&& c):a(std::move(c)){}
	free_algebra(const free_algebra<R>& p) :a(p.a) {}
	free_algebra(int c) :a(1,R(c)) {}

	static free_algebra _0()
	{
		return free_algebra();
	}
	static free_algebra _1()
	{
		return free_algebra(1);
	}
	int degree() const
	{
		return a.size() - 1;
	}

	free_algebra operator-() const
	{
		free_algebra p(*this);
		for (auto& s : p.a)
			s = -s;
		return p;
	}

	template<typename H=R>
	H operator()(const H& u)
	{
		H r = 0,w=1;
		for (const auto& s : a)
		{
			r += w*s;
			w *= u;
		}
		return r;
	}

	const free_algebra& operator+() const
	{
		return *this;
	}

	free_algebra& operator+=(const free_algebra &p)
	{
		for (int i = 0; i <= degree(); i++)
			if (i > p.degree())
				break;
			else a[i] += p.a[i];
		for (int i = degree() + 1; i <= p.degree(); i++)
			a.push_back(p.a[i]);
		reduce();
		return *this;
	}
	free_algebra& operator+=(const R& p)
	{
		return *this+=free_algebra(p);
	}
	free_algebra& operator-=(const R& p)
	{
		return *this += free_algebra(p);
	}

	free_algebra& operator*=(const R& p)
	{
		std::for_each(a.begin(), a.end(), [&p](auto& v) {v *= p; });
		reduce();
		return *this;
	}
	free_algebra& operator-=(const free_algebra& p)
	{
		for (int i = 0; i <= degree(); i++)
			if (i > p.degree())
				break;
			else a[i] -= p.a[i];
		for (int i = degree() + 1; i <= p.degree(); i++)
			a.push_back(p.a[i]);
		reduce();
		return *this;
	}
	free_algebra & operator*=(const free_algebra& p)
	{
		if ((degree() < 0) || (p.degree() < 0))
		{
			*this = _0();
			return *this;
		}
		int m = degree() + p.degree();
		free_algebra q;
		q.a.resize(m + 1);
		for (int i = 0; i <= degree(); i++)
			for (int j = 0; j <= p.degree(); j++)
				q.a[i + j] += a[i] * p.a[j];
		*this =(q);
		reduce();
		return *this;
	}
	const R& coeff(int n) const
	{
		if (n > degree())
			return R::_0();
		return a.at(n);
	}
	virtual ring& operator+=(int n) 
	{
		return *this += R(n);
	}
	virtual ring& operator-=(int n)
	{
		return *this -= R(n);
	}
	virtual ring& operator*=(int n)
	{
		return *this *= R(n);
	}
	bool is_zero() const
	{
		return a.empty();
	}
	bool is_one() const
	{
		return !a.empty() && a.at(0).is_one();
	}
protected:
	void reduce() {
		while (!a.empty()&&(a.back() == R::_0()))
			a.pop_back();
	}
	std::vector<R> a;
};

template<typename R>
free_algebra<R> operator+(const free_algebra<R>& a, const free_algebra<R>& b)
{
	free_algebra<R> p(a);
	return p += b;
}


template<typename R>
free_algebra<R> operator-(const free_algebra<R>& a, const free_algebra<R>& b)
{
	free_algebra<R> p(a);
	return p -= b;
}

template<typename R>
free_algebra<R> operator*(const free_algebra<R>& a, const free_algebra<R>& b)
{
	free_algebra<R> p(a);
	return p *= b;
}

template<typename R>
free_algebra<R> operator+(const free_algebra<R>& a, const R& b)
{
	free_algebra<R> p(a);
	return p += b;
}


template<typename R>
free_algebra<R> operator-(const free_algebra<R>& a, const R& b)
{
	free_algebra<R> p(a);
	return p -= b;
}

template<typename R>
free_algebra<R> operator*(const free_algebra<R>& a, const R& b)
{
	free_algebra<R> p(a);
	return p *= b;
}

template<typename R>
free_algebra<R> operator+(const R& b,const free_algebra<R>& a)
{
	free_algebra<R> p(a);
	return p += b;
}


template<typename R>
free_algebra<R> operator-(const R& b,const free_algebra<R>& a)
{
	free_algebra<R> p(a);
	return p -= b;
}

template<typename R>
free_algebra<R> operator*(const R& b,const free_algebra<R>& a)
{
	free_algebra<R> p(a);
	return p *= b;
}

template<typename R>
std::ostream& operator<<(std::ostream& H, const free_algebra<R>& p)
{
	H << "(";
	for (int i = 0; i <= p.degree(); i++)
		if (i < p.degree())
			H << p.coeff(i) << ", ";
		else H << p.coeff(i) <<")";
	return H;
}