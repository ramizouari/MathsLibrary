#pragma once
#include "polynomial.h"
namespace math_rz
{
	template<typename F>
	polynomial<F> lagrange_basis_vector(const std::vector<F>& x,int k)
	{
		int m = x.size();
		polynomial<F> e_k = 1;
		for (int i = 0; i < m; i++)
			if (i == k)
				continue;
			else e_k *= polynomial<F>({ -x[i],1 }) / F(x[k] - x[i]);
		return e_k;
	}

	template<typename F>
	polynomial<F> lagrange(const std::vector<F>& x, const std::vector<F>& y)
	{
		if (x.size() != y.size())
			throw std::exception("x & y must have the same size");
		int m = x.size();
		polynomial<F> L = 0;
		for (int i = 0; i < m; i++)
			L += y[i] * lagrange_basis_vector(x, i);
		return L;
	}
}