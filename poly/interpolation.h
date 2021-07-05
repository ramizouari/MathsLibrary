#pragma once
#include "polynomial.h"
#include "linalg/transformation/fft/cooley_tuckey.h"
namespace math_rz::poly
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

	template<typename K>
	polynomial<K> newton_interpolation(const std::vector<K>& x, const std::vector<K>& y)
	{
		if (x.size() != y.size())
			throw std::exception("x & y must have the same size");
		std::vector<std::vector<K>> divided_difference(x.size());
		divided_difference[0] = y;
		int n = y.size() - 1;
		for (int r = 1; r <= n; r++)
		{
			divided_difference[r].resize(n-r+1);
			for (int i = 0; i + r <= n; i++)
				divided_difference[r][i] = (divided_difference[r - 1][i + 1] - divided_difference[r - 1][i]) / (x[i + r] - x[i]);
		}
		polynomial<K> L=0;
		for (int r = n; r >= 0; r--)
			L = L * polynomial<K>({ -x[r],1 }) + divided_difference[r][0];
		return L;
	}

	template<unsigned int n>
	polynomial<complex> fft_interpolation(const linalg::finite_dimensional_vector_space<complex, n>& y)
	{
		static linalg::fft::dynamic_cooley_tuckey<true> CT(n);
		auto R = CT(y.get_vect());
		for(auto &x:R)
			x /= (int)n;
		return (R);
	}

}