#pragma once
#include "integer.h"
#include <type_traits>
#include "field.h"

namespace math_rz
{
	template<typename long long n,bool is_prime=false>
	class cyclic:public std::conditional<is_prime,field,ring>::type
	{
	public:
		cyclic(integer s=0):w((s+n)%n){}
		cyclic(long long k):w((k+n)%n){}
		cyclic& operator+=(const cyclic& b)
		{
			w += b.w;
			w%=n;
			return *this;
		}
		cyclic& operator-=(const cyclic& a)
		{
			w += (n - a.w);
			w %= n;
			return *this;
		}
		cyclic& operator*=(const cyclic& a)
		{
			w *= a.w;
			w %= n;
			return *this;
		}


		cyclic& operator+=(int b)
		{
			return *this+=cyclic(b);
		}
		cyclic& operator-=(int a)
		{
			return *this-=cyclic(a);
		}
		cyclic& operator*=(int a)
		{
			return *this*=cyclic(a);
		}
		cyclic operator-() const
		{
			return cyclic(n - w);
		}
		operator integer& () { return w; }
		operator const integer& () const{ return w; }
		operator long long()const { return w._v; }
		bool is_zero() const
		{
			return w.is_zero();
		}
		bool is_one() const
		{
			return w.is_one();
		}
		cyclic inv()const
		{
			return math_rz::bezout(w,n).first;
		}

		cyclic operator/=(const cyclic& a)
		{
			return *this *= a.inv();
		}
	private:
		integer w;
	};
	template<long long n,bool is_prime>
	cyclic<n, is_prime> operator+(const cyclic<n, is_prime>& a, const cyclic<n, is_prime>& b)
	{
		auto c(a);
		return c += b;
	}
	template<long long n,bool is_prime >
		cyclic<n,is_prime> operator-(const cyclic<n, is_prime>& a, const cyclic<n, is_prime>& b)
	{
		auto c(a);
		return c -= b;
	}
	template<long long n,bool is_prime>
	cyclic<n, is_prime> operator*(const cyclic<n, is_prime>& a, const cyclic<n, is_prime>& b)
	{
		auto c(a);
		return c *= b;
	}

	template<long long n, bool is_prime>
	cyclic<n, is_prime> operator/(const cyclic<n, is_prime>& a, const cyclic<n, is_prime>& b)
	{
		auto c(a);
		return c /= b;
	}
	template<long long n>
	using cyclic_ring= cyclic<n, false>;
	template<long long n>
	using cyclic_field = cyclic<n, true>;
}

