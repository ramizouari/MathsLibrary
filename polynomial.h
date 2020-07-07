#pragma once
#include "free_algebra.h"
#include "integral_ring.h"
#include "rational_extension.h"
template<typename F>
class polynomial:virtual public free_algebra<F>,virtual public integral_ring
{
public:
	using free_algebra<F>::a;
	polynomial() {};
	polynomial(const free_algebra<F> &p):free_algebra<F>(p){}
	polynomial(free_algebra<F>&& p) :free_algebra<F>(std::move(p)) {}
	polynomial(const std::vector<F>& p) :free_algebra<F>(p) {}
	polynomial(std::vector<F>&& p) :free_algebra<F>(std::move(p)) {}
	polynomial(const F &p):free_algebra<F>(p){}
	polynomial(int c):free_algebra<F>(c){}
	bool operator!=(const polynomial& p) const
	{
		if (p.degree() != this->degree())
			return true;
		for (int i = 0; i < p.degree(); i++)
			if (p.a.at(i) != this->a.at(i))
				return true;
		return false;
	}
	bool operator==(const polynomial& p) const
	{
		return !(*this != p);
	}
	static polynomial _0()
	{
		return free_algebra<F>::_0();
	}
	static polynomial _1()
	{
		return free_algebra<F>::_1();
	}
	static std::pair<polynomial, polynomial> euclidean_division(const polynomial& p, const polynomial& q) 
	{
		if (p.degree() < q.degree())
			return std::make_pair(F::_0(), p);
		polynomial r(p);
		int m(r.degree()), n(q.degree());
		polynomial s;
		s.a.resize(m - n + 1);
		for (; m >= n;m--)
		{
			F k(r.a.at(m) / q.a.at(n));
			s.a.at(m-n) = k;
			if (k == F::_0())
			{
				r.a.pop_back();
				continue;
			}
			for (int i = 1; i <= n; i++)
				r.a.at(m - i) -= k * q.a.at(n - i);
			r.a.pop_back();
		}
		r.reduce();
		return std::make_pair(s, r);
	}
	static polynomial gcd(const polynomial& p, const polynomial& q)
	{
		if (p.degree() < q.degree())
			return gcd(q, p);
		std::pair<polynomial,polynomial> R(euclidean_division(p, q));
		if (R.second.a.empty())
			return q.normalize();
		else return gcd(q, R.second);
	}
	polynomial div(const polynomial& q) const
	{
		return euclidean_division(*this, q).first;
	}
	polynomial mod(const polynomial& q) const
	{
		return euclidean_division(*this, q).second;
	}
	polynomial& operator+=(const polynomial& p)
	{
		this->free_algebra<F>::operator+=(p);
		return *this;
	}
	polynomial& operator-=(const polynomial& p)
	{
		this->free_algebra<F>::operator-=(p);
		return *this;
	}
	polynomial& operator*=(const polynomial& p)
	{
		this->free_algebra<F>::operator*=(p);
		return *this;
	}
	polynomial& operator+=(const F& p)
	{
		return this->free_algebra<F>::operator+=(p);
	}
	polynomial& operator-=(const F& p)
	{
		return this->free_algebra<F>::operator-=(p);
	}
	polynomial& operator*=(const F& p)
	{
		return this->free_algebra<F>::operator*=(p);
	}
	polynomial& operator/=(const F& p)
	{
		std::for_each(a.begin(), a.end(), [&p](auto& v) {v /= p; });
		this->reduce();
		return *this;
	}
	polynomial normalize() const
	{
		polynomial p(*this);
		auto dominant_coeff = p.a.at(p.degree());
		p /= dominant_coeff;
		return p;
	}
};

template<typename F>
polynomial<F> operator+(const polynomial<F>& a, const polynomial<F>& b)
{
	polynomial<F> p(a);
	return p += b;
}


template<typename F>
polynomial<F> operator-(const polynomial<F>& a, const polynomial<F>& b)
{
	polynomial<F> p(a);
	return p -= b;
}

template<typename F>
polynomial<F> operator*(const polynomial<F>& a, const polynomial<F>& b)
{
	polynomial<F> p(a);
	return p *= b;
}

template<typename F>
polynomial<F> operator+(const polynomial<F>& a, const F& b)
{
	polynomial<F> p(a);
	return p += b;
}


template<typename F>
polynomial<F> operator-(const polynomial<F>& a, const F& b)
{
	polynomial<F> p(a);
	return p -= b;
}

template<typename F>
polynomial<F> operator*(const polynomial<F>& a, const F& b)
{
	polynomial<F> p(a);
	return p *= b;
}

template<typename F>
polynomial<F> operator+(const F& b, const polynomial<F>& a)
{
	polynomial<F> p(a);
	return p += b;
}


template<typename F>
polynomial<F> operator-(const F& b, const polynomial<F>& a)
{
	polynomial<F> p(a);
	return p -= b;
}

template<typename F>
polynomial<F> operator*(const F& b, const polynomial<F>& a)
{
	polynomial<F> p(a);
	return p *= b;
}

template <typename F>
using rational_function = rational_extension<polynomial<F>>;