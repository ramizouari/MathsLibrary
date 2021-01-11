#pragma once
#include "norm.h"
namespace math_rz::poly::structure
{
	template<typename F>
	class inner_product_topology:public norm_topology<F>
	{
	public:
		virtual F inner_product(const polynomial<F>& p, const polynomial<F>& q) const = 0;
		virtual real_field norm(const polynomial<F>& p) const
		{
			return std::sqrt(static_cast<real_field>(inner_product(p, p)));
		}
	};

	template<typename F>
	class L2_vect_inner_product:public inner_product_topology<F>
	{
	public:
		F inner_product(const polynomial<F>& p, const polynomial<F>& q) const
		{
			if (p.degree() <= q.degree())
			{
				auto w = p.conj();
				return std::inner_product(w.get_vect().cbegin(), w.get_vect().cend(),
					q.get_vect().cbegin(), F(0));
			}
			else
			{
				auto w = q.conj();
				return std::inner_product(q.get_vect().cbegin(), q.get_vect().cend(),
					p.get_vect().cbegin(), F(0)).conj();
			}

		}
	};


	template<typename F>
	class L2_function_inner_product :public inner_product_topology<F>
	{
		std::shared_ptr<integrator<F, F>> I_ptr;
	public:
		L2_function_inner_product(std::shared_ptr<integrator<F, F>> _I_ptr) :I_ptr(_I_ptr) {}
		L2_function_inner_product(integrator<F, F>* _I_ptr) :I_ptr(_I_ptr) {}
		F inner_product(const polynomial<F>& p, const polynomial<F>& q) const
		{
			return I_ptr->integrate
			(
				general_function<F, F>([&](const F& u)->F
					{
						return p(u).conj()*q(u);
					})
			);
		}
	};
}