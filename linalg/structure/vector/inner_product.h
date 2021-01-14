#pragma once
#include "norm.h"
namespace math_rz::linalg::structure::vector
{
	template<typename K,int n>
	class inner_product_topology:public norm_topology<K,n>
	{
		using vector_type = norm_topology<K, n>::vector_type;
	public:
		virtual K inner_product(const vector_type& p, const vector_type& q) const = 0;
		virtual K dot_product(const vector_type& p, const vector_type& q) const
		{
			return inner_product(p.conj(), q);
		}
		virtual real_field norm(const vector_type& p) const
		{
			return std::sqrt(static_cast<real_field>(inner_product(p, p)));
		}
	};

	template<typename K,int n>
	class L2_vect_inner_product:public inner_product_topology<K,n>
	{
		using vector_type = norm_topology<K, n>::vector_type;
	public:
		K inner_product(const vector_type& p, const vector_type& q) const
		{
			auto w = p.conj();
			return std::inner_product(w.get_vect().cbegin(), w.get_vect().cend(),
				q.get_vect().cbegin(), K(0));
		}
	};

	template<typename K, int n>
	class L2_induced_vect_inner_product :public inner_product_topology<K, n>
	{
		using vector_type = norm_topology<K, n>::vector_type;
		square_matrix<K, n> M;
	public:
		L2_induced_vect_inner_product(const square_matrix<K, n>& P):M(P){}
		K inner_product(const vector_type& p, const vector_type& q) const
		{
			auto w = p.conj(),z=(M*q);
			return std::inner_product(w.get_vect().cbegin(), w.get_vect().cend(),
				z.get_vect().cbegin(), K(0));
		}
	};
}