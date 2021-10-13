#pragma once
#include "characterestic_polynomial.h"
#include "poly/interpolation.h"
namespace math_rz::linalg::characterestic
{
	template<ring_constraints::field R, int n>
	class interpolation_method :public characterestic_polynomial<R, n>
	{
	protected:
		using polynomial_type = std::conditional_t<ring_constraints::field<R>, poly::polynomial<R>, poly::free_algebra<R>>;
	public:
		virtual polynomial_type operator()(const matrix<R, n, n>& A) const
		{
			matrix<R,n,n> M = A;
			std::vector<R> X(n+1), Y(n+1);
			for (int i = 0; i <= n; i++)
			{
				X[i] = i;
				Y[i] = M.det();
				for (int j = 0; j < n; j++)
					M[j][j] = M[j][j] - 1;
			}
			return poly::newton_interpolation(X, Y);
		}
	};

}