#pragma once
#include "linalg/matrix.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "function.h"

namespace math_rz {
	template<typename F, int n, int m, int p = 2>
	class derivator
	{
		using E1 = math_rz::Lp_finite_dimensional_space<F, p, m>;
		using E2 = math_rz::Lp_finite_dimensional_space<F, p, n>;
	public:
		derivator(E1 p0,F _eps) :x0(p0),eps(_eps) {}

		matrix<F, n, m> jacobian(const function<E1,E2> & f) const
		{
			E1 s = x0;
			matrix<F, n, m> M;
			for (int i = 0; i < n; i++)
			{
				s[i] += eps;
				F k = F(1) / eps;
				E2 h = k*(f(s) - f(x0));
				s[i] -= eps;
				for (int j = 0; j < m; j++)
					M[j][i] = h[j];

				
			}
			return M;
		}

	private:
		E1 x0;
		F eps;
	};
}