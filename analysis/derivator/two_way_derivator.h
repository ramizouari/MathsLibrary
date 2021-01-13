#pragma once
#include "derivator.h"

namespace math_rz::analysis
{

	/*
	* General Case:
	*/
	template<typename E1,
		typename E2>
		requires linalg::vector_space_constraint::vector_space_over_same_base_field<E1, E2>
		class two_way_derivator : public derivator<E1, E2>
	{
	private:
		using F = typename derivator<E1, E2>::base_field;
		static inline constexpr int n = derivator<E1, E2>::codomain_dimension;
		static inline constexpr int m = derivator<E1, E2>::domain_dimension;
	public:
		two_way_derivator(F _eps) :eps(_eps) {}

		derivator<E1, E2>::matrix_type jacobian(const function<E1, E2>& f, const E1& x0) const override
		{
			E1 s1 = x0,s2=x0;
			typename derivator<E1, E2>::matrix_type M;
			F k = F(.5) / eps;
			if constexpr (derivator<E1, E2>::domain_dimension > 1 &&
				derivator<E1, E2>::codomain_dimension > 1)
			{
				for (int i = 0; i < m; i++)
				{
					s1[i] += eps;
					s2[i] -= eps;
					E2 h = k * (f(s1) - f(s2));
					s1[i] -= eps;
					s2[i] += eps;
					for (int j = 0; j < n; j++)
						M[j][i] = h[j];

				}
			}
			else if constexpr (derivator<E1, E2>::codomain_dimension > 1)
			{
				s1 += eps;
				s2 -= eps;
				E2 h = k * (f(s1) - f(s2));
				s1 -= eps;
				s2 += eps;
				for (int j = 0; j < n; j++)
					M[j][0] = h[j];
			}
			else if constexpr (derivator<E1, E2>::domain_dimension > 1)
			{
				for (int i = 0; i < m; i++)
				{
					s1[i] += eps;
					s2[i] -= eps;
					E2 h = k * (f(s1) - f(s2));
					s1[i] -= eps;
					s2[i] += eps;
					M[0][i] = h;
				}
			}
			else
			{
				s1 += eps;
				s2 -= eps;
				E2 h = k * (f(s1) - f(s2));
				s1 -= eps;
				s2 += eps;
				M[0][0] = h;
			}
			return M;
		}

	private:
		F eps;
	};


}