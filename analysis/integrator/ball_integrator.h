#pragma once
#include "integrator.h"
#include <numbers>
#include "complex.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/general_function.h"
#include "special_integrator.h"
namespace math_rz::analysis {

	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F> requires (E::dimension == 3)
	class ball_integrator :
		public  special_integrator<E, F, E, F>
	{
	public:
		ball_integrator(integrator<E, F>* _I)
			:special_integrator<E, F, E, F>(_I)
		{

		}

		ball_integrator(std::shared_ptr<integrator<E, F>> _I)
			:special_integrator<E, F, E, F>(_I)
		{

		}

		F integrate(const function<E, F>& f) const override
		{
			return this->I_ptr->integrate
			(
				general_function<E, F>([&](const E& s)->F 
					{
						E u({ s[0] * std::cos(s[1]) * std::sin(s[2]),
							s[0] * std::sin(s[1]) * std::sin(s[2]),
							s[0] *std::cos(s[2]) });
						return std::pow(s[0], 2) * std::sin(s[2]) * f(u);
					})
			);
		}
	};
}