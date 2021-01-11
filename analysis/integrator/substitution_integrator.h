#pragma once
#include "integrator.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "analysis/function.h"
#include "special_integrator.h"
#include "analysis/derivator/derivator.h"

namespace math_rz
{
	/*
	* This class is used for a change of variable when integrating multivariable function
	*/
	template<typename E,typename F>
	class substitution_integrator :public  special_integrator<E, F, E, F>
	{
	
		function<E, E>& phi;
		std::shared_ptr<derivator<E, E>> D_ptr;
	public:
		substitution_integrator(function<E, E>& _mu,integrator<E, F> *_I_ptr,
			derivator<E,E> *_D_ptr)
			:phi(_mu),special_integrator<E, F,E, F>(_I_ptr),D_ptr(_D_ptr)
		{

		}

		substitution_integrator(function<E, E>& _mu, std::shared_ptr<integrator<E, F>> _I_ptr,
			derivator<E, E>* _D_ptr)
			:phi(_mu), special_integrator<E, F, E, F>(_I_ptr), D_ptr(_D_ptr)
		{

		}

		substitution_integrator(function<E, E>& _mu, integrator<E, F>* _I_ptr,
			std::shared_ptr<derivator<E, E>> _D_ptr)
			:phi(_mu), special_integrator<E, F, E, F>(_I_ptr), D_ptr(_D_ptr)
		{

		}

		substitution_integrator(function<E, E>& _mu, std::shared_ptr<integrator<E, F>> _I_ptr,
			std::shared_ptr<derivator<E, E>> _D_ptr)
			:phi(_mu), special_integrator<E, F, E, F>(_I_ptr), D_ptr(_D_ptr)
		{

		}

		F integrate(const function<E, F>& f) const override
		{
			return this->I_ptr->integrate(general_function<E, F>([&](const E& u)
				{
					auto K = D_ptr->jacobian_det(phi,u).abs();
					return K * f(phi(u));
				}
			));
		}
	};
}