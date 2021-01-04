#pragma once
#include <complex>
#include "absalg/field.h"
#include "real_field.h"
namespace math_rz
{
	class complex :public std::complex<real_field>, public field
	{
	public:
		complex(const std::complex<real_field>& a) :std::complex<real_field>(a) {}
		complex(const std::complex<double>& a) :std::complex<real_field>(a) {}
		complex(const std::complex<float>& a) :std::complex<real_field>(a) {}
		complex(const std::complex<long double>& a) :std::complex<real_field>(a) {}

		complex(const real_field& a = 0, const real_field& b = 0);
		complex(double s);
		complex(float s);
		complex(long double s);
		complex(int a);
		static complex _0() { return std::complex<real_field >(); }
		static complex _1() { return std::complex<real_field>(1); }
		bool is_zero() const;
		bool is_one() const;
		real_field abs() const;
		complex conj() const;
		inline operator std::complex<real_field>() { return *this; }
		complex& operator/=(const real_field& s);
		complex operator-() const;
		complex& operator*=(int a);
		complex& operator*=(const std::complex<real_field> &s);
		//inline operator real_field() { return real(); }
		complex& operator/=(int a);
	};

	complex operator""_c(long double a);
	complex operator""_c(unsigned long long a);
	complex operator*(const complex& s, int a);
	complex operator*(int a,const complex& s);

	//complex operator/(const complex &s ,const real_field& a);

}

