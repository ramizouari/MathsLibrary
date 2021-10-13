#pragma once
#include "characterestic_polynomial.h"
#include "poly/interpolation.h"
namespace math_rz::linalg::characterestic
{
	template<ring_constraints::ring R, int n>
	class faddev_leverrier :public characterestic_polynomial<R, n>
	{
	protected:
		using polynomial_type = std::conditional_t<ring_constraints::field<R>, poly::polynomial<R>, poly::free_algebra<R>>;
	public:
		virtual polynomial_type operator()(const matrix<R, n, n>& A) const
		{
			std::vector<R> S(n + 1);
			S[n] = 1;
			matrix<R, n, n> C;
			for (int i = n - 1; i >= 0; i--)
			{
				for (int j = 0; j < n; j++)
					C[j][j] += S[i + 1];
				C = A * C;
				S[i] = -C.tr() / R(n - i);
			}
			return poly::polynomial<R>(S);
		}
	};
}