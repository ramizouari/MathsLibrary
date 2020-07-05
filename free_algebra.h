#pragma once
#include "ring.h"
#include <vector>
#include <algorithm>
#include <iostream>

template<typename R>
class free_algebra: virtual public ring
{
public:
	free_algebra(){}
	free_algebra(R m):a(1,m){}
	free_algebra(std::vector<R>&& c):a(std::move(c)){}
	free_algebra(const free_algebra<R>& p) :a(p.a) {}

	static free_algebra _0()
	{
		return free_algebra();
	}
	static free_algebra _1()
	{
		return free_algebra(1);
	}
	const free_algebra& I0() const
	{
		return _0();
	}
	const free_algebra& I1() const
	{
		return _1();
	}
	int degree() const
	{
		return a.size() - 1;
	}
	free_algebra& operator+=(const free_algebra &p)
	{
		for (int i = 0; i <= degree(); i++)
			if (i > p.degree())
				break;
			else a[i] += p.a[i];
		for (int i = degree() + 1; i <= p.degree(); i++)
			a.push_back(p.a[i]);
		sync();
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
		sync();
		return *this;
	}
	free_algebra& operator*=(const free_algebra& p)
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
		*this = std::move(q);
		sync();
		return *this;
	}
	R coeff(int n) const
	{
		if (n > degree())
			return R::_0();
		return a.at(n);
	}
protected:
	void sync() {
		while (!a.empty()&&(a.front() == R::_0()))
			a.pop_back();
	}
private:
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
std::ostream& operator<<(std::ostream& H, const free_algebra<R>& p)
{
	H << "(";
	for (int i = 0; i <= p.degree(); i++)
		if (i < p.degree())
			H << p.coeff(i) << ", ";
		else H << p.coeff(i) <<")";
	return H;
}