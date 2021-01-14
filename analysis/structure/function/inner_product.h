#pragma once
#include "norm.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::analysis::structure::function
{
	template<typename A, math_rz::linalg::vector_space_constraint::vector_space B>
	class inner_product_topology:public norm_topology<A,B>
	{
	public:
		using K = typename norm_topology<A,B>::K;
		using function_type = math_rz::analysis::function<A,B>;
	public:
		virtual K inner_product(const function_type& p, const function_type& q) const = 0;
		virtual K dot_product(const function_type& p, const function_type& q) const
		{
			return inner_product(p.conj(), q);
		}

		virtual real_field norm(const function_type& p) const
		{
			return std::sqrt(static_cast<real_field>(inner_product(p, p)));
		}
	};

	template<typename A,typename B>
	class L2_function_inner_product:public inner_product_topology<A,B>
	{
		using function_type = math_rz::analysis::function<A,B>;
		using K = typename  norm_topology<A, B>::K;
		std::shared_ptr<integrator<A, K>> I_ptr;
	public:
		L2_function_inner_product(integrator<A, K>* _I_ptr) : I_ptr(_I_ptr) {}
		L2_function_inner_product(std::shared_ptr<integrator<A, K>> _I_ptr): I_ptr(_I_ptr) {}		
		
		K inner_product(const function_type& p, const function_type& q) const
		{
			return I_ptr->integrate
			(
				general_function<A, K>([&](const K& u)->K
					{
						return p(u).inner_product(q(u));
					})
			);
		}
	};

	/*template<typename A,typename B>
	class L2_induced_vect_inner_product :public inner_product_topology<A,B>
	{
		using function_type = norm_topology<A,B>::function_type;
		square_matrix<A,B> M;
	public:
		L2_induced_vect_inner_product(const square_matrix<A,B>& P):M(P){}
		F inner_product(const function_type& p, const function_type& q) const
		{
			auto w = p.conj(),z=(M*q);
			return std::inner_product(w.get_vect().cbegin(), w.get_vect().cend(),
				z.get_vect().cbegin(), F(0));
		}
	};*/
}