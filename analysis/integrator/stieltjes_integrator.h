#pragma once
#include "integrator.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "analysis/function.h"

namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F>
	class stieltjes_integrator :public  integrator<real_field, F,real_field>
	{
		integer cuts;
		real_field a, b;
		function<real_field, E> &mu;
	public:
		stieltjes_integrator(function<real_field,E>&_mu,real_field _a, real_field _b, integer _cuts = 100) 
			:mu(_mu),a(_a), b(_b), cuts(_cuts)
		{

		}
		real_field integrate(const function<real_field, F>& f) const override
		{
			real_field eps = (b - a) / cuts;
			real_field u = a,v=u+eps, result;
			for (int i = 0; i < cuts; i++, u += eps, v += eps)
				if constexpr (F::dimension == 1 && E::dimension == 1)
					result += (mu(v) - mu(u)) * f(u);
				else if constexpr (std::is_same<E, F>::value)
					 result += E(mu(v) - mu(u)).inner_product(f(u));
				else result += E(mu(v) - mu(u)).norm()*f(u);
			return result;
		}
	};
}