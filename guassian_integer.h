#pragma once
#include <complex>
#include "integer.h"
#include "absalg/integral_ring.h"

namespace math_rz
{
	class guassian_integer :public std::complex<integer>, public ring
	{
	public:
		guassian_integer(const std::complex<integer>& a) :std::complex<integer>(a) {}
		guassian_integer(const integer& a = 0, const integer& b = 0);
		guassian_integer(int a);
		guassian_integer(const guassian_integer& a, int b);
		static guassian_integer _0() { return std::complex<integer>(); }
		static guassian_integer _1() { return std::complex<integer>(1); }
		bool is_zero() const;
		bool is_one() const;
		inline operator std::complex<integer>() { return *this; }
		guassian_integer conj() const;
	};

}