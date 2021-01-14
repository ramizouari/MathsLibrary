#pragma once
#pragma once
#include "integrator.h"
#include <numbers>
#include "complex.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/general_function.h"
#include "special_integrator.h"
namespace math_rz::analysis {
	/*
	* This class is responsible for integrating a function over a sphere of a given radius R (equals 1 by default)
	*/
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F> requires (E::dimension == 3)
	class spherical_integrator :
		public  special_integrator<E, F,
		linalg::coordinate_space<typename E::base_field,2>,F>
	{
		using E2 = linalg::coordinate_space<typename E::base_field, 2>;
		real_field R;
	public:
		spherical_integrator(integrator<E2, F>* _I, real_field r = 1) 
			:special_integrator<E,F,E2,F>(_I), R(r)
		{

		}

		spherical_integrator(std::shared_ptr<integrator<E2, F>> _I, real_field r = 1) 
			:special_integrator<E, F, E2, F>(_I), R(r)
		{

		}

		F integrate(const function<E, F>& f) const override
		{
			return this->I_ptr->integrate
			(
				general_function<E2, F>([&](const E2& s)->F 
					{
						E u({ R*std::cos(s[0]) * std::sin(s[1]),
							R*std::sin(s[0]) * std::sin(s[1]),
							R*std::cos(s[1]) });
						return std::pow(R,2)*std::sin(s[1])*f(u);
					})
			);
		}
	};
}