#pragma once
#pragma once
#include "integrator.h"
#include <numbers>
#include "complex.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/general_function.h"
#include "special_integrator.h"
namespace math_rz {
	/*
	* This class is responsible for integrating a function over a sphere of a given radius R (equals 1 by default)
	*/
	template<typename F ,int p=2>
	class spherical_integrator :
		public  special_integrator<Lp_finite_dimensional_space<real_field,p,3>, F,
		Lp_finite_dimensional_space<real_field, p, 2>,F>
	{
		using Lp = Lp_finite_dimensional_space<real_field, p, 3>;
		using Lp2 = Lp_finite_dimensional_space<real_field, p, 2>;
		real_field R;
	public:
		spherical_integrator(integrator<Lp2, F>* _I, real_field r = 1) 
			:special_integrator<Lp,F,Lp2,F>(_I), R(r)
		{

		}

		spherical_integrator(std::shared_ptr<integrator<Lp2, F>> _I, real_field r = 1) 
			:special_integrator<Lp, F, Lp2, F>(_I), R(r)
		{

		}

		F integrate(const function<Lp, F>& f) const override
		{
			return this->I_ptr->integrate
			(
				general_function<Lp2, F>([&](const Lp2& s)->F 
					{
						Lp u({ R*std::cos(s[0]) * std::sin(s[1]),
							R*std::sin(s[0]) * std::sin(s[1]),
							R*std::cos(s[1]) });
						return std::pow(R,2)*std::sin(s[1])*f(u);
					})
			);
		}
	};
}