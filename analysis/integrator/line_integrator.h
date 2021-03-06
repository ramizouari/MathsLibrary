#pragma once
#include "special_integrator.h"
#include "analysis/general_function.h"
#include "analysis/derivator/derivator.h"
#include "stieltjes_integrator.h"

namespace math_rz::analysis
{
	/*
	* This class is responsible for calculating the integral of a function along a (smooth enough) curve
	* There are 2 cases:
	* 1. The function is scalar
	* 2. The function is a vector function
	*/
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F> 
	requires (E::dimension == F::dimension || F::dimension == 1)
	class line_integrator : public special_integrator<E, F, real_field, real_field,real_field>
	{
		const function<real_field, E>& phi_ptr;
		std::shared_ptr<derivator<real_field, E>> D;
	public:
		line_integrator(const function<real_field,E>&_phi_ptr,integrator<real_field,real_field,real_field>*_I_ptr,
			derivator<real_field, E>*_D)
			:special_integrator<E,F,real_field,real_field, real_field>(_I_ptr),phi_ptr(_phi_ptr),D(_D){}

		line_integrator(const function<real_field, E>& _phi_ptr, 
			std::shared_ptr<integrator<real_field, real_field, real_field>> _I_ptr,
			derivator<real_field, E>* _D)
			:special_integrator<E, F, real_field, real_field, real_field>(_I_ptr), phi_ptr(std::move(_phi_ptr)), D(_D) {}

		line_integrator(const function<real_field, E>& _phi_ptr, integrator<real_field, real_field, real_field>* _I_ptr,
			std::shared_ptr<derivator<real_field, E>> _D)
			:special_integrator<E, F, real_field, real_field, real_field>(_I_ptr), phi_ptr(_phi_ptr), D(_D) {}

		line_integrator(const function<real_field, E>& _phi_ptr, 
			std::shared_ptr<integrator<real_field, real_field, real_field>> _I_ptr,
			std::shared_ptr<derivator<real_field, E>> _D)
			:special_integrator<E, F, real_field, real_field, real_field>(_I_ptr), phi_ptr(std::move(_phi_ptr)), D(_D) {}
		
		real_field integrate(const function<E, F>& f) const override
		{
			if constexpr(F::dimension == 1)
				return this->I_ptr->integrate
				(
					general_function<real_field, F>([&](const real_field& s)->F 
						{
							return D->derivative(phi_ptr,s).norm()*f((phi_ptr)(s));
						})
				);
			else return this->I_ptr->integrate
				(
					general_function<real_field, real_field>([&](const real_field& s)->real_field 
						{
							return f((phi_ptr)(s)).inner_product(D->derivative(phi_ptr, s))  ;
						})
				);
		}
	};

	
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F>
	class stieltjes_line_integrator : public special_integrator<E, F, real_field, F, real_field>
	{
		function<real_field, E>& phi_ptr;
		std::shared_ptr<stieltjes_integrator<E,F>> I_ptr;
	public:
		stieltjes_line_integrator(function<real_field, E>& _phi_ptr,
			std::shared_ptr<stieltjes_integrator<E,F>> _I_ptr)
			:I_ptr(_I_ptr),
			special_integrator<E, F, real_field, F, real_field>(I_ptr), phi_ptr(_phi_ptr) 
		{}
		stieltjes_line_integrator(function<real_field, E>& _phi_ptr,
			stieltjes_integrator<E,F>* _I_ptr)
			:I_ptr(_I_ptr),
			special_integrator<E, F, real_field, F, real_field>(I_ptr), phi_ptr(_phi_ptr)
		{}

		real_field integrate(const function<E, F>& f) const override
		{
			if constexpr (E::dimension == 1)
				return this->I_ptr->integrate
				(
					general_function<real_field, F>([&](const real_field& s)->F 
						{
							return f((phi_ptr)(s));
						})
				);
			else return this->I_ptr->integrate
			(
				general_function<real_field, F>([&](const real_field& s)->F
					{
						return f((phi_ptr)(s));
					})
			);
		}
	};
}