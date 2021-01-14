#pragma once
#include "special_integrator.h"
#include "analysis/general_function.h"
#include "analysis/derivator/derivator.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "linalg/multiplicator/cross_product.h"
#include "prob/uniform_real_generator.h"
#include "linalg/eigen.h"

namespace math_rz::analysis
{
	/*
	* This integrator gives the integral along the boundary of a (smooth enough) manifold parameterized by phi
	*/
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F>
	requires (F::dimension == E::dimension || F::dimension == 1)
		class boundary_integrator : public special_integrator<E, 
		F,
		std::conditional_t<E::dimension==2,typename E::base_field,
		linalg::coordinate_space<typename E::base_field, E::dimension-1>>, 
		real_field, 
		real_field>
	{
		using base_field = typename E::base_field;
		inline static constexpr int n = E::dimension;
		using H = std::conditional_t<n==2,base_field,linalg::coordinate_space<base_field, n-1>>;
		function<H, E>& phi_ptr;
		std::shared_ptr<derivator<H, E>> D;
		bool reverse_orientation;
	public:
		boundary_integrator(function<H, E>& _phi_ptr, integrator<H, real_field, real_field>* _I_ptr,
			derivator<H, E>* _D, bool reverse = false)
			:special_integrator<E, F, H, real_field, real_field>(_I_ptr), phi_ptr(_phi_ptr), D(_D),
			reverse_orientation(reverse)
		{}

		boundary_integrator(function<H, E>& _phi_ptr, std::shared_ptr<integrator<H, real_field, real_field>> _I_ptr,
			std::shared_ptr<derivator<H, E>> _D, bool reverse = false)
			:special_integrator<E, F, H, real_field, real_field>(_I_ptr), phi_ptr(_phi_ptr), D(_D),
			reverse_orientation(reverse) {}

		boundary_integrator(function<H, E>&& _phi_ptr, integrator<H, real_field, real_field>* _I_ptr,
			derivator<H, E>* _D)
			:special_integrator<E, F, H, real_field, real_field>(_I_ptr),
			phi_ptr(std::move(_phi_ptr)), D(_D) {}

		real_field integrate(const function<E, F>& f) const override
		{
			if constexpr (F::dimension == 1)
				return this->I_ptr->integrate
				(
					general_function<H, F>([&](const H& s)->F
						{
							auto differential = D->jacobian(phi_ptr, s);
							/*
							* The determinant of the following matrix gives the  squared product of the singular values of the differential
							*/
							auto K = std::sqrt(square_matrix<base_field, n-1>(differential.H() * differential).det());
							return  K * f(phi_ptr(s));
						})
				);
			else return this->I_ptr->integrate
			(
				general_function<H, real_field>([&](const H& s)->real_field
					{
						/*
						* Get the transpose of the differential
						*/
						auto differential_T = D->jacobian(phi_ptr, s).T();
						/*
						* Store the differential in p
						*/
						std::vector<std::vector<base_field>> p=differential_T.get_vect_vect();
						static uniform_real_generator G(-1,1,432);
						/*
						* Add to p a random vector
						* Almost surely, the random won't be in the span of p before the addition
						*/
						p.push_back(G.generate_vector<n>().get_vect());
						/*
						* construct a matrix from p
						*/
						linalg::matrix<base_field, n, n> M(std::move(p));
						/*
						* if M is singular, return 0
						* Almost surely, M is singular iff the differential is not injective
						*/
						if (M.nullity() > 0)
							return 0;
						/*
						* Else apply the Gram-Schmidt process to the matrix
						* The last vector w in this process will be orthogonal to the hyperplane defined by the
						* differential, or the hyperplane defined by the differential is tangent to the manifold
						*-> The vector w will be a unit normal vector to the manifold
						*/
						linalg::gram_schmidt_inplace<base_field,n,n>(M);
						/*
						* To get a consistent choice of the unit normal vector, we will impose that the resultant
						* gram-schmidt basis must be positively oriented
						*/
						bool positive_orientation = square_matrix<base_field, n>(M).det() > 0;
						/*
						* the vector w which is the unit normal vector to the manifold
						*/
						linalg::coordinate_space<base_field, n> w(std::move(M[n - 1]));
						if (reverse_orientation ^ positive_orientation)
							w = -w;
						/*
						* Calculate the differential change of hypervolume
						*/
						auto K = std::sqrt(square_matrix<base_field,n-1>(differential_T * differential_T.H()).det());
						return  K*w.inner_product(f((phi_ptr)(s)));
					})
			);
		}
	};
}