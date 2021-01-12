#pragma once
#include "polynomial.h"
namespace math_rz
{
	template<typename K>
	polynomial<K> lagrange_basis_vector(const std::vector<K>& x,int k)
	{
		int m = x.size();
		polynomial<K> e_k = 1;
		for (int i = 0; i < m; i++)
			if (i == k)
				continue;
			else e_k *= polynomial<K>({ -x[i],1 }) / K(x[k] - x[i]);
		return e_k;
	}

	template<typename K>
	polynomial<K> lagrange(const std::vector<K>& x, const std::vector<K>& y)
	{
		if (x.size() != y.size())
			throw std::exception("x & y must have the same size");
		int m = x.size();
		polynomial<K> L = 0;
		for (int i = 0; i < m; i++)
			L += y[i] * lagrange_basis_vector(x, i);
		return L;
	}
}