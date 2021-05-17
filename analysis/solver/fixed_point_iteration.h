#pragma once
#include "analysis/function.h"
#include "real_field.h"
#include "fixed_point_finder.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::metric_space E>
	class fixed_point_iteration:public fixed_point_finder<E>
	{
		E x0;
		real_field eps=1e-5;
	public:
		fixed_point_iteration(const E& _x0) :x0(_x0) {}

		virtual E fixed_point(const function<E, E>& f) const override
		{
			E x = x0;
			while (f(x).distance(x) > eps) x = f(x);
			return x;
		}
	};
}