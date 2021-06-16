#pragma once
#include "constrained_minimiser.h"
namespace math_rz::analysis
{
	template<typename E,typename F>
	class lagrangian_minimiser:public constrained_minimiser<E>
	{
		using H = linalg::vector_space_constraint::product_space<E, F>;
		minimiser<H>& M;
		function<E, F>& C;
	public:
		lagrangian_minimiser(function<E,F>&_C,minimiser<H> &_M):C(_C),M(_M){}
		virtual E argmin(const function<E, real_field>& f) const
		{
			general_function<H, real_field> L([&](const H& u)->real_field
				{
					E x=u.get<0,E::dimension>();
					F µ=u.get<E::dimension,H::dimension>();
					return f(x) + µ.inner_product(C(x));
				});
			auto v = M.argmin(L);
			return v.get<0, E::dimension>();
		}

	
	};
}