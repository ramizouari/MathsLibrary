#pragma once
#pragma once
#include "special_integrator.h"
#include <numbers>
#include "complex.h"

namespace math_rz::analysis {

	template<typename F >
	class complex_circular_integrator :public  special_integrator<complex, F,real_field,F>
	{
		real_field R;
	public:
		complex_circular_integrator(std::shared_ptr<integrator<real_field,F>> _I, real_field r=1) 
			:special_integrator<complex, F, real_field, F>(_I),R(r)
		{

		}

		complex_circular_integrator(integrator<real_field, F>* _I, real_field r = 1)
			:special_integrator<complex, F, real_field, F>(_I), R(r)
		{

		}

		F integrate(const function<complex, F>& f) const override
		{
			return this->I_ptr->integrate
			(
				general_function<real_field,F>([&](const real_field& s)->F 
					{
						real_field p = s;
						auto w = complex(R*std::cos(p), R * std::sin(p));
						return complex(0,R)* w * f(w);
					})
			);
		}
	};

	template< linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F > requires (E::dimension ==2)
	class circular_integrator 
		:public  special_integrator<E, F,real_field,F>
	{
		real_field R;
	public:
		circular_integrator(std::shared_ptr<integrator<real_field, F>> _I, real_field r = 1) 
			:special_integrator<E, F, real_field, F>(_I), R(r)
		{

		}

		F integrate(const function<E, F>& f) const override
		{
			return this->I_ptr->integrate
			(
				general_function<real_field, F>([&](const real_field& s)->F 
					{
						auto w = E({ R * std::cos(s), R * std::sin(s) });
						return R * f(w);
					})
			);
		}
	};
}