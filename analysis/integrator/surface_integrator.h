#pragma once
#include "special_integrator.h"
#include "analysis/general_function.h"
#include "analysis/derivator/derivator.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "linalg/multiplicator/cross_product.h"

namespace math_rz::analysis
{
	/*
	* This class is responsible for integrating a function f along a (smooth enough) surface parameterized by phi
	* There are two cases for the integral
	* 1. The function is scalar
	* 
	* 2. The function is a vector field (on the euclidean space of dimension 3), 
	* for this case the surface must be orientable so that the notion of integral makes sense
	*/
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F> 
	requires (E::dimension == F::dimension && E::dimension == 3 || F::dimension == 1)
	class surface_integrator : public special_integrator<E, F, 
		linalg::coordinate_space<typename E::base_field,2>, real_field, real_field>
	{
		using E2 = linalg::coordinate_space<typename E::base_field, 2>;
		function<E2, E>& phi_ptr;
		std::shared_ptr<derivator<E2, E>> D;
		bool reverse_orientation;
	public:
		surface_integrator(function<E2, E>& _phi_ptr, integrator<E2, real_field, real_field>* _I_ptr,
			derivator<E2, E>* _D,bool reverse=false)
			:special_integrator<E, F, E2, real_field, real_field>(_I_ptr), phi_ptr(_phi_ptr), D(_D),
			reverse_orientation(reverse)
		{}

		surface_integrator(function<E2, E>& _phi_ptr, std::shared_ptr<integrator<E2, real_field, real_field>> _I_ptr,
			std::shared_ptr<derivator<E2, E>> _D, bool reverse = false)
			:special_integrator<E, F, E2, real_field, real_field>(_I_ptr), phi_ptr(_phi_ptr), D(_D),
		reverse_orientation(reverse){}

		surface_integrator(function<E2, E>&& _phi_ptr, integrator<E2, real_field, real_field>* _I_ptr,
			derivator<E2, E>* _D)
			:special_integrator<E, F, E2, real_field, real_field>(_I_ptr), 
			phi_ptr(std::move(_phi_ptr)), D(_D) {}
		
		real_field integrate(const function<E, F>& f) const override
		{
			if constexpr (F::dimension == 1)
				return this->I_ptr->integrate
				(
					general_function<E2, F>([&](const E2& s)->F 
					{
						//Used to calculate cross product
						static linalg::cross_product<F> C;
						E u, v;
						auto differential = D->jacobian(phi_ptr, s);
						/*
						* The determinant of the following matrix gives the product of the singular values of differential squared
						*/
						auto K = std::sqrt(square_matrix<typename E2::base_field,2>(differential.H() * differential).det());
						return  K * f(phi_ptr(s));
					})
				);
			else return this->I_ptr->integrate
			(
				general_function<E2, real_field>([&](const E2& s)->real_field 
				{
					//Used to calculate cross product
					static linalg::cross_product<typename E::base_field> C;
					E u, v;
					auto differential = D->jacobian(phi_ptr, s);
					for (int i = 0; i < u.dimension; i++)
					{
						u[i] = differential[i][0];
						v[i] = differential[i][1];
					}
					auto w = C.multiply(u,v);
					if (reverse_orientation)
						w = -w;
					return  w.inner_product(f((phi_ptr)(s)));
				})
			);
		}
	};
}