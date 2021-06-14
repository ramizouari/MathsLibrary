#pragma once
#include "analysis/function.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "analysis/derivator/derivator.h"
namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::normed_vector_space E>
	class halley_method
	{
		derivator<E, E>& D;
		E x0;
		real_field eps = 1e-5;
	public:
		halley_method(E _x0, derivator<E, E>& d) :D(d), x0(_x0) {}
		virtual E root(const function<E, E>& f) const
		{
			E x = x0;
			while (f(x).norm() > eps)
				if constexpr (E::dimension > 1)
					x = x - D.jacobian(f, x).inv() * f(x);
				else x = x - f(x) / D.derivative(f, x);
			return x;
		}
	};
}