#pragma once
#include "gradient_descent.h"
namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::inner_product_space E ,
		linalg::vector_space_constraint::inner_product_space F>
	class constrained_barzilai_borwein_gradient_descent :public gradient_descent<E>
	{
		mutable E s;
		function<E, F>& I;
	public:
		constrained_barzilai_borwein_gradient_descent(const E& _x0, derivator<E, real_field>& d,
			function<E,F> &_I, real_field _p)
			:gradient_descent<E>(_x0, d),I(_I)
		{
			this->p = _p;
		}

		E argmin(const function<E, real_field>& f) const override
		{
			this->p = 0.1;
			s = this->x0;
			E x = s - this->p * this->D.gradient(f, s);
			for (; (x-s).norm() > this->eps; x -= this->p * this->D.gradient(f, x))
			{
				update_rate(f, x);
				bool feasible = true;
				s = x;
			}
			return x;
		}

		virtual void update_rate(const function<E, real_field>& f, const E& x)const override
		{
			auto L = this->D.gradient(f, x) - this->D.gradient(f, s);
			this->p = (L).inner_product(x - s) / L.inner_product(L);
		}
	};

	template<typename E>
	using BB_gradient_descent = barzilai_borwein_gradient_descent<E>;
}