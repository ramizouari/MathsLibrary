#pragma once
#include "analysis/function.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "analysis/derivator/derivator.h"
#include "linalg/square_matrix.h"
namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::inner_product_space E>
	class secant_method
	{
		E x0,x1;
		real_field eps = 1e-5;
		derivator<E, E>& D;
	public:
		secant_method(E _x0,E _x1,derivator<E,E> &d) : D(d),x0(_x0),x1(_x1) {}
		virtual E root(const function<E, E>& f) const
		{
			E p = x0,q=x1,tmp;
			linalg::square_matrix<typename E::base_field,E::dimension> invJ=D.jacobian(f,x0).inv();
			while (f(p).norm() > eps)
			{
				E delta_x = (q - p),delta_f=(f(q)-f(p));
				invJ = invJ + (delta_x - invJ * delta_f).outer_product( delta_x) * invJ/((delta_f).inner_product(invJ*delta_x));
				tmp = q;
				q = q - invJ*f(q);
				p = tmp;
			}
			return q;
		}
	};

	template<linalg::vector_space_constraint::inner_product_space E> requires (E::dimension == 1)
	class secant_method<E>
	{
		E x0, x1;
		real_field eps = 1e-5;
	public:
		secant_method(E _x0, E _x1) : x0(_x0), x1(_x1) {}
		virtual E root(const function<E, E>& f) const
		{
			E p = x0, q = x1, tmp;
			while (f(p).norm() > eps)
			{
				tmp = q;
				q = q - f(q) * (q - p) / (f(q) - f(p));
				p = tmp;
			}
			return q;
		}
	};


}