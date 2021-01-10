#pragma once
#include "special_integrator.h"
#include "analysis/general_function.h"
#include "analysis/derivator/derivator.h"

namespace math_rz
{
	template<typename E,typename F>
	class line_integrator : public special_integrator<E, F, real_field, real_field,real_field>
	{
		function<real_field, E>& phi_ptr;
		std::shared_ptr<derivator<real_field, 1, E::dimension, real_field, E>> D;
	public:
		line_integrator(function<real_field,E>&_phi_ptr,integrator<real_field,real_field,real_field>*_I_ptr,
			derivator<real_field, 1, E::dimension, real_field, E>*_D)
			:special_integrator<E,F,real_field,real_field, real_field>(_I_ptr),phi_ptr(_phi_ptr),D(_D){}

		line_integrator(function<real_field, E>&& _phi_ptr, integrator<real_field, real_field, real_field>* _I_ptr,
			derivator<real_field, 1, E::dimension, real_field, E>* _D)
			:special_integrator<E, F, real_field, real_field, real_field>(_I_ptr), phi_ptr(std::move(_phi_ptr)), D(_D) {}
		real_field integrate(const function<E, F>& f) const override
		{
			if constexpr(F::dimension == 1)
			return this->I_ptr->integrate(general_function<real_field, F>([&](const real_field& s)->F {
				return D->derivative(phi_ptr,s).norm()*f((phi_ptr)(s));
				}));
			else return this->I_ptr->integrate(general_function<real_field, real_field>([&](const real_field& s)->real_field {
				return f((phi_ptr)(s)).inner_product(D->derivative(phi_ptr, s))  ;
				}));
		}
	};
}