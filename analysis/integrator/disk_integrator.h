#pragma once
#include "integrator.h"
#include <numbers>
#include "complex.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/general_function.h"
#include "special_integrator.h"
namespace math_rz::analysis {

	/*
	* This class is responsible for integrating a function over (possibly portion of) a disk
	* If the integral is over the whole disk, the integrator must be compatible with the standard mesure
	* over the set [0,R]x[0,2pi] where R is the radius of the disk
	*/
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F>
	requires (E::dimension == 2)
	class disk_integrator :
		public  special_integrator<E, F,E, F>
	{
		
	public:
		disk_integrator(integrator<E, F>* _I, real_field r = 1)
			:special_integrator<E, F, E, F>(_I)
		{

		}

		disk_integrator(std::shared_ptr<integrator<E, F>> _I, real_field r = 1)
			:special_integrator<E, F, E, F>(_I)
		{

		}

		F integrate(const function<E, F>& f) const override
		{
			return this->I_ptr->integrate
			(
				general_function<E, F>([&](const E& s)->F 
					{
						E u({ s[0] * std::cos(s[1]), s[0] * std::sin(s[1])});
						return s[0] * f(u);
					})
			);
		}
	};
}