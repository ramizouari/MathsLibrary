#pragma once
#include "gradient_descent.h"
namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::normed_vector_space E>
	class fixed_rate_gradient_descent :public gradient_descent<E>
	{
	public:
		fixed_rate_gradient_descent(const E& _x0, derivator<E, real_field>& d,real_field _p) 
			:gradient_descent<E>(_x0,d) 
		{
			this->p = _p;
		}

		virtual void update_rate(const function<E, real_field>& f, const E& x)const override {}
	};
}