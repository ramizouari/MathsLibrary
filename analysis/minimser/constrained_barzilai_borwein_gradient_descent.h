#pragma once
#include "gradient_descent.h"
namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::inner_product_space E ,
		linalg::vector_space_constraint::inner_product_space F,
		linalg::vector_space_constraint::inner_product_space H>
	class constrained_barzilai_borwein_gradient_descent :public gradient_descent<E>
	{
		mutable E s;
		function<F, H>& I;
		bool constraint_feasible(const E& x) const
		{
			bool feasible = true;
			I(x.get<0,F::dimension>()).foreach([&feasible](const auto& a)
				{
					if (a > 0)
						feasible = false;
				});
			return feasible;
		}
	public:
		constrained_barzilai_borwein_gradient_descent(const E& _x0, derivator<E, real_field>& d,
			function<F,H> &_I, real_field _p)
			:gradient_descent<E>(_x0, d),I(_I)
		{
			this->p = _p;
		}

		E argmin(const function<E, real_field>& f) const override
		{
			this->p = 0.1;
			s = this->x0;
			E grad = this->D.gradient(f, s);
			E x = s - this->p * grad;
			bool feasible;
			for (; (x-s).norm() > this->eps; grad=this->D.gradient(f,s),x -= this->p * grad)
			{
				while (!constraint_feasible(x) && this->p>this->eps)
				{
					this->p /= 2;
					x = s - this->p * grad;
				}
				update_rate(f, x);
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

	template<typename E,typename F,typename H>
	using CBB_gradient_descent = constrained_barzilai_borwein_gradient_descent<E,F,H> ;
}