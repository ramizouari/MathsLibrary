#pragma once
#include <complex>
#include "integer.h"
class guassian_integer:public std::complex<integer>
{
public:
	guassian_integer();
	guassian_integer(const integer& a, const integer& b);
};

