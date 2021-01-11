#pragma once
#include "integrator.h"
#include <numbers>
#include "complex.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/general_function.h"
#include "special_integrator.h"
namespace math_rz {

	template<typename F, int p = 2>
	class disk_integrator :
		public  special_integrator<Lp_finite_dimensional_space<real_field, p, 2>, F,
		Lp_finite_dimensional_space<real_field, p, 2>, F>
	{
		using Lp = Lp_finite_dimensional_space<real_field, p, 2>;
		using Lp2 = Lp;
	public:
		disk_integrator(integrator<Lp2, F>* _I, real_field r = 1)
			:special_integrator<Lp, F, Lp2, F>(_I)
		{

		}

		disk_integrator(std::shared_ptr<integrator<Lp2, F>> _I, real_field r = 1)
			:special_integrator<Lp, F, Lp2, F>(_I)
		{

		}

		F integrate(const function<Lp, F>& f) const override
		{
			return this->I_ptr->integrate
			(
				general_function<Lp2, F>([&](const Lp2& s)->F 
					{
						Lp u({ s[0] * std::cos(s[1]), s[0] * std::sin(s[1])});
						return s[0] * f(u);
					})
			);
		}
	};
}