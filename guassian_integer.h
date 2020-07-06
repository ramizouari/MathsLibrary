#pragma once
#include <complex>
#include "integer.h"
#include "integral_ring.h"
class guassian_integer :public std::complex<integer>, public ring
{
public:
	guassian_integer(const std::complex<integer>& a) :std::complex<integer>(a) {}
	guassian_integer(const integer& a=0, const integer& b=0);
	static guassian_integer _0() { return std::complex<integer>(); }
	static guassian_integer _1() { return std::complex<integer>(1); }
	operator std::complex<integer>() { return *this; }
	
};

