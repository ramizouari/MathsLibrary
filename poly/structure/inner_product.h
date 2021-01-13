#pragma once
#include "norm.h"
namespace math_rz::poly::structure
{
	template<typename K>
	class inner_product_topology:public norm_topology<K>
	{
	public:
		virtual K inner_product(const polynomial<K>& p, const polynomial<K>& q) const = 0;
		virtual real_field norm(const polynomial<K>& p) const
		{
			return std::sqrt(static_cast<real_field>(inner_product(p, p)));
		}
	};

	template<typename K>
	class L2_vect_inner_product:public inner_product_topology<K>
	{
	public:
		K inner_product(const polynomial<K>& p, const polynomial<K>& q) const
		{
			if (p.degree() <= q.degree())
			{
				auto w = p.conj();
				return std::inner_product(w.get_vect().cbegin(), w.get_vect().cend(),
					q.get_vect().cbegin(), K(0));
			}
			else
			{
				auto w = q.conj();
				return std::inner_product(q.get_vect().cbegin(), q.get_vect().cend(),
					p.get_vect().cbegin(), K(0)).conj();
			}

		}
	};


	template<typename K>
	class L2_function_inner_product :public inner_product_topology<K>
	{
		std::shared_ptr<analysis::integrator<K, K>> I_ptr;
	public:
		L2_function_inner_product(std::shared_ptr<analysis::integrator<K, K>> _I_ptr) :I_ptr(_I_ptr) {}
		L2_function_inner_product(analysis::integrator<K, K>* _I_ptr) :I_ptr(_I_ptr) {}
		K inner_product(const polynomial<K>& p, const polynomial<K>& q) const
		{
			return I_ptr->integrate
			(
				general_function<K, K>([&](const K& u)->K
					{
						return p(u).conj()*q(u);
					})
			);
		}
	};
}