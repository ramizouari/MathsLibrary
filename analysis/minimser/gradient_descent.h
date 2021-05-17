#pragma once
#include "minimiser.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "analysis/derivator/derivator.h"

namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::normed_vector_space E>
	class gradient_descent : public minimiser<E>
	{
	protected:
		E x0;
		mutable real_field p=.1;
		real_field eps=1e-5;
		derivator<E, real_field>& D;
	public:
		gradient_descent(const E& _x0,derivator<E,real_field> &d):x0(_x0),D(d) {}
		void set_seed(const E& _x0) { x0 = _x0; }
		virtual E argmin(const function<E, real_field>& f) const override
		{
			E x = x0;
			for (; D.gradient(f, x).norm() > eps; x -= p * D.gradient(f, x))
				update_rate(f, x);
			return x;
		}

		virtual void update_rate(const function<E, real_field>& f, const E& x)const = 0;
	};
}