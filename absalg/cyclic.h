#pragma once
#include "integer.h"
#include <type_traits>
#include "field.h"
#include <optional>

namespace math_rz
{
    class dynamic_cyclic
    {
    public:
        explicit dynamic_cyclic(integer _n, integer s = 0);
        dynamic_cyclic(const dynamic_cyclic &o,integer s):dynamic_cyclic(o.n,s){}
        dynamic_cyclic& operator+=(const dynamic_cyclic& b);
        dynamic_cyclic& operator-=(const dynamic_cyclic& a);
        dynamic_cyclic& operator*=(const dynamic_cyclic& a);
        dynamic_cyclic& operator+=(int b);
        dynamic_cyclic& operator-=(int a);
        dynamic_cyclic& operator*=(int a);
        dynamic_cyclic operator-() const;
        operator integer& ();
        operator const integer& () const;
        operator long long()const;
        bool is_zero() const;
        bool is_one() const;
        dynamic_cyclic inv()const;
        dynamic_cyclic operator/=(const dynamic_cyclic& a);
        static dynamic_cyclic primitive_root(int m);
        static dynamic_cyclic primitive_unity_root(int n,int m);
        static dynamic_cyclic principal_unity_root(int n, int m);
    private:
        integer n;
        integer w;
    };
    dynamic_cyclic operator+(const dynamic_cyclic& a, const dynamic_cyclic& b);
    dynamic_cyclic operator-(const dynamic_cyclic& a, const dynamic_cyclic& b);
    dynamic_cyclic operator*(const dynamic_cyclic& a, const dynamic_cyclic& b);
    dynamic_cyclic operator/(const dynamic_cyclic& a, const dynamic_cyclic& b);

    template<long long n,bool is_prime=false>
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
			if constexpr (is_prime)
				return math_rz::pow(*this, n - 2);
			else return math_rz::bezout<integer>(w, n).first;
		}

		cyclic operator/=(const cyclic& a)
		{
			return *this *= a.inv();
		}
		inline static cyclic primitive_root()
		{
			static std::optional<cyclic> root;
			if(!root.has_value())
				root.emplace((std::int64_t)dynamic_cyclic::primitive_root(n));
			return root.value();
		}
		inline static cyclic primitive_unity_root(int m)
		{
			return (std::int64_t)dynamic_cyclic::primitive_unity_root(n, m);
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
	namespace finite_fields
	{
		using F2 = cyclic_field<2>;
		using F3 = cyclic_field<3>;
		using F5 = cyclic_field<5>;
		using F7 = cyclic_field<7>;
		using F11 = cyclic_field<11>;
		using F13 = cyclic_field<13>;
		using F17 = cyclic_field<17>;
	}
	namespace finite_rings
	{
		using namespace finite_fields;
		using Z4 = cyclic_ring<4>;
		using Z8 = cyclic_ring<8>;
		using Z6 = cyclic_ring<6>;
		using Z9 = cyclic_ring<9>;
		using Z15 = cyclic_ring<15>;
	}
	
}

