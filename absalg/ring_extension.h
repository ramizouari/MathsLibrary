#pragma once
#include <integer.h>
#include <real_field.h>
#include <type_traits>
#include "poly/polynomial.h"
#include "field.h"
namespace math_rz 
{
	template <long long x>
	concept is_quadratic_extension = x==2;
	template<typename R,bool is_field,long long ...k>
	class ring_extension :public std::conditional<is_field, field, ring>::type
	{
	public:
		ring_extension(const std::vector<R>& a) :p(a) { reduce(); }
		ring_extension(int s = 0) :p(s) { reduce(); }
		const static inline polynomial<R> extension_polynomial 
			= polynomial<R>(std::vector<R>{ k... });
		ring_extension(const R& a) :p(a) { reduce(); }
		ring_extension& operator+=(const ring_extension &w)
		{
			p += w.p;
			return *this;
		}

		ring_extension& operator-=(const ring_extension& w)
		{
			p -= w.p;
			return *this;
		}
		ring_extension& operator*=(const ring_extension& w)
		{
			p *= w.p;
			reduce();
			return *this;
		}
		ring_extension& operator/=(const ring_extension& w)
		{
			p *= w.p.inv();
			return *this;
		}

		ring_extension inv() const
		{
			return math_rz::bezout<polynomial<ring_extension>>(p, extension_polynomial).first;
		}
		bool is_zero() const override
		{
			return p.is_zero();
		}
		bool is_one() const override
		{
			return p.is_one();
		}
		const polynomial<R>& get_polynomial() const { return p; }
		inline static constexpr long long extension_degree = sizeof ...(k) - 1;
		ring_extension quad_conj() const requires is_quadratic_extension<extension_degree>
		{
			if (p.degree() < 2)
				return *this;
			return ring_extension({ p.coeff(0),-p.coeff(1) });
		}
	private:
		void reduce()
		{
			if (p.degree() >= extension_polynomial.degree())
				p = p.mod(extension_polynomial);
		}
		polynomial<R> p;
	};

	template<typename R, bool is_field,long long ...k>
	ring_extension<R,is_field,k...> operator+(const ring_extension<R,is_field,k...>& p, 
		const ring_extension<R,is_field,k...>& q)
	{
		ring_extension<R,is_field,k...>h = p;
		return h += q;
	}

	template<typename R, bool is_field,long long ...k>
	ring_extension<R,is_field,k...> operator-(const ring_extension<R, is_field,k...>& p,
		const ring_extension<R, is_field,k...>& q)
	{
		ring_extension<R,is_field,k...>h = p;
		return h -= q;
	}

	template<typename R, bool is_field,long long ...k>
	ring_extension< R,is_field, k... > operator*(const ring_extension<R,is_field,k...>& p, 
		const ring_extension<R,is_field,k...>& q)
	{
		ring_extension<R,is_field,k...>h = p;
		return h *= q;
	}

	template<typename R, bool is_field,long long ...k>
	ring_extension<R,is_field,k...> operator/(const ring_extension<R, is_field,k...>& p, 
		const ring_extension<R,is_field,k...>& q)
	{
		ring_extension<R,is_field,k...>h = p;
		return h /= q;
	}

	template<typename R, bool is_field,long long ...k>
	std::ostream&  operator<<(std::ostream&H, const ring_extension<R,is_field,k...>& q)
	{
		return H << q.get_polynomial();
	}

	using eisenstein_rational = ring_extension<rational_extension<integer>,true,1, 1, 1>;
	using eisenstein_real = ring_extension<real_field,true, 1, 1, 1>;
	using eisenstein_integer = ring_extension<integer,false, 1, 1, 1>;
	
	template<long long n>
	using quadratic_field = ring_extension<rational_extension<integer>,n<0, -n,0, 1>;

	template<long long n>
	using cubic_field = ring_extension<rational_extension<integer>,false, -n, 0, 0, 1>;

}