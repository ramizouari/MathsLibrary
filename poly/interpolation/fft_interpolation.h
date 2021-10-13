#pragma once
#include "linalg/transformation/fft/cooley_tuckey.h"
#include "poly/polynomial.h"
#include "complex.h"

namespace poly::interpolation
{

	template<int n>
	polynomial<complex> fft_interpolation(const linalg::finite_dimensional_vector_space<complex, n>& y)
	{
		static linalg::fft::cooley_tuckey<std::bit_ceil<n>> CT;
		return CT(CT(CT(y))).get_vect();
	}

    template<int p,bool prime,int n>
    polynomial<complex> fft_interpolation(const linalg::finite_dimensional_vector_space<cyclic<p,prime>, n>& y)
    {
        static linalg::fft::cooley_tuckey<std::bit_ceil<n>> CT;
        return CT(CT(CT(y))).get_vect();
    }

}