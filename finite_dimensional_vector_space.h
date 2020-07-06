#pragma once
#include "vector_space.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>

template<typename F,int n>
class finite_dimensional_vector_space:public vector_space<F>
{
public:
	finite_dimensional_vector_space(const std::vector<F>& a):u(n)
	{
		if (n != a.size())
			throw std::domain_error("Dimensions are not compatible");
		std::copy(a.begin(), a.end(), u.begin());
	}
	finite_dimensional_vector_space(std::vector<F>&& a) :u(std::move(a))
	{
		if (n != u.size())
			throw std::domain_error("Dimensions are not compatible");
	}
	finite_dimensional_vector_space():u(n){}
	finite_dimensional_vector_space& operator+=(const finite_dimensional_vector_space& o)
	{
		for (int i = 0; i < n; i++)
			u.at(i) += o.u.at(i);
		return *this;
	}
	finite_dimensional_vector_space& operator-=(const finite_dimensional_vector_space& o)
	{
		for (int i = 0; i < n; i++)
			u.at(i) -= o.u.at(i);
		return *this;
	}
	finite_dimensional_vector_space& operator*=(const F& k)
	{
		for (int i = 0; i < n; i++)
			u.at(i) *= k;
		return *this;
	}

	finite_dimensional_vector_space& operator/=(const F& k)
	{
		for (int i = 0; i < n; i++)
			u.at(i) /= k;
		return *this;
	}

	finite_dimensional_vector_space& operator*=(int k)
	{
		for (int i = 0; i < n; i++)
			u.at(i) *=k;
		return *this;
	}
	finite_dimensional_vector_space& operator/=(int k)
	{
		for (int i = 0; i < n; i++)
			u.at(i) /= k;
		return *this;
	}

	finite_dimensional_vector_space operator-() const
	{
		finite_dimensional_vector_space p;
		std::transform(u.begin(), u.end(), p.u.begin(), [](auto a) {return -a; });
		return p;
	}
	const F& operator[](int i) const
	{
		return u[i];
	}
	F& operator[](int i)
	{
		return u[i];
	}
	const F& at(int i) const
	{
		return u.at(i);
	}
	F& at(int i)
	{
		return u.at(i);
	}
protected:
	std::vector<F> u;
};
template <typename F, int n>
using coordinate_space = finite_dimensional_vector_space<F, n>;

template <typename F,int n>
finite_dimensional_vector_space<F, n> operator+(
	const finite_dimensional_vector_space<F, n>& a,const finite_dimensional_vector_space<F, n>& b)
{
	auto c(a);
	return c += b;
}

template <typename F, int n>
finite_dimensional_vector_space<F, n> operator-(
	const finite_dimensional_vector_space<F, n>& a, const finite_dimensional_vector_space<F, n>& b)
{
	auto c(a);
	return c -= b;
}

template <typename F, int n>
finite_dimensional_vector_space<F, n> operator*(
	const F& k, const finite_dimensional_vector_space<F, n>& a)
{
	auto c(a);
	return c *= k;
}

template <typename F, int n>
finite_dimensional_vector_space<F, n> operator/(
	const F& k, const finite_dimensional_vector_space<F, n>& a)
{
	auto c(a);
	return c /= k;
}

template <typename F, int n>
finite_dimensional_vector_space<F,n> operator*(int k, const finite_dimensional_vector_space<F, n>& a)
{
	auto c(a);
	return c *= k;
}

template <typename F, int n>
finite_dimensional_vector_space<F, n> operator/(int k, const finite_dimensional_vector_space<F, n>& a)
{
	auto c(a);
	return c /= k;
}

template <typename F,int n>
std::ostream& operator<<(std::ostream& H, const finite_dimensional_vector_space<F, n>& p)
{
	H << "( ";
	for (int i = 0; i < n; i++)
		if (i == n-1) H << p.at(i) << " )";
		else H << p.at(i) << ", ";
	return H;
}