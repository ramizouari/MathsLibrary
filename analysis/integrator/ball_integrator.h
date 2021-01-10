#pragma once
#include "integrator.h"
#include <numbers>
#include "complex.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/general_function.h"
#include "special_integrator.h"
namespace math_rz {

	template<typename F, int p = 2>
	class ball_integrator :
		public  special_integrator<Lp_finite_dimensional_space<real_field, p, 3>, F,
		Lp_finite_dimensional_space<real_field, p, 3>, F>
	{
		using Lp = Lp_finite_dimensional_space<real_field, p, 3>;
		using Lp3 = Lp;
	public:
		ball_integrator(integrator<Lp3, F>* _I)
			:special_integrator<Lp, F, Lp3, F>(_I)
		{

		}

		ball_integrator(std::shared_ptr<integrator<Lp3, F>> _I)
			:special_integrator<Lp, F, Lp3, F>(_I)
		{

		}

		F integrate(const function<Lp, F>& f) const override
		{
			return this->I_ptr->integrate(general_function<Lp3, F>([&](const Lp3& s)->F {
				Lp u({ s[0] * std::cos(s[1]) * std::sin(s[2]),
					s[0] * std::sin(s[1]) * std::sin(s[2]),
					s[0] *std::cos(s[2]) });
				return std::pow(s[0], 2) * std::sin(s[2]) * f(u);
				}));
		}
	};
}