#pragma once
#include <complex>
#include "absalg/field.h"
#include "real_field.h"
namespace math_rz
{
	class complex :public std::complex<real_field>, public field
	{
	public:
		using base_field = complex;
		inline constexpr static int dimension = 1;
		complex(const std::complex<real_field>& a) :std::complex<real_field>(a) {}
		complex(const std::complex<double>& a) :std::complex<real_field>(a) {}
		complex(const std::complex<float>& a) :std::complex<real_field>(a) {}
		complex(const std::complex<long double>& a) :std::complex<real_field>(a) {}
		inline complex(integer s):std::complex<real_field>(s._v) {};
		complex(const real_field& a = 0, const real_field& b = 0);
		complex(double s);
		complex(float s);
		complex(long double s);
		complex(int a);
		complex(long long a);
		static complex _0() { return std::complex<real_field >(); }
		static complex _1() { return std::complex<real_field>(1); }
		bool is_zero() const;
		bool is_one() const;
		real_field abs() const;
		complex conj() const;
		complex inner_product(const complex&z)const;
		complex dot_product(const complex& z)const;
		complex inv()const;
		real_field norm() const;
		inline operator std::complex<long double>() { return std::complex<long double>(real(),imag()); }
		//inline operator std::complex<double>() { return std::complex<double>(real(), imag()); }
		explicit inline operator real_field() { return real(); }
		complex& operator/=(const real_field& s);
		complex operator-() const;
		complex& operator*=(int a);
		//complex& operator*=(long double a);
		//complex& operator*=(double a);
		//complex& operator*=(real_field);
		complex& operator*=(const std::complex<real_field>& s);
		complex& operator/=(const std::complex<real_field>& s);
		//inline operator real_field() { return real(); }
	};

	complex operator""_c(long double a);
	complex operator""_c(unsigned long long a);
//	complex operator*(const complex& s, int a);
	//complex operator*(const complex& s, const complex& a);
//	complex operator*(int a,const complex& s);
	//complex operator*(double a, const complex& s);
	//complex operator*(long double a, const complex& s);
	//complex operator*(const real_field &a, const complex& b);
//	complex operator/(const std::complex<real_field>& a, const std::complex<real_field>& b);
	//complex operator/(const complex &s ,const real_field& a);
}

