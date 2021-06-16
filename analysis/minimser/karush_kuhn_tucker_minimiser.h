#pragma once
#include "constrained_minimiser.h"
#include "absalg/direct_product.h"
namespace math_rz::analysis
{
	template<typename E, typename F,typename G>
	class karush_kuhn_tucker_minimiser :public constrained_minimiser<E>
	{
		using H = linalg::vector_space_constraint::product_space<linalg::vector_space_constraint::product_space<E, F>,G>;
		minimiser<H>& M;
		function<E, F>& C;
		function<E, G>& I;

	public:
		karush_kuhn_tucker_minimiser(function<E, F>& _C, function<E, G>& _I, minimiser<H>& _M) :C(_C), M(_M),I(_I) {}
		virtual E argmin(const function<E, real_field>& f) const override
		{
			general_function<H, real_field> L([&](const H& u)->real_field
				{
					E x = u.get<0, E::dimension>();
					F µ = u.get<E::dimension, E::dimension+F::dimension>();
					G λ = u.get<E::dimension + F::dimension, H::dimension>();
					λ.foreach([](auto& a) {a *= a; });
					return f(x)+ µ.inner_product(C(x))+ λ.inner_product(I(x));
				});
			auto v = M.argmin(L);
			return v.get<0, E::dimension>();
		}


	};

	template<typename E,typename F,typename G>
	using KKT_minimiser = karush_kuhn_tucker_minimiser<E, F, G>;
}