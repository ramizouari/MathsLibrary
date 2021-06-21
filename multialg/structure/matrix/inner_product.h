#pragma once
#include "norm.h"
namespace math_rz::linalg::structure::matrix
{
	template<typename K,int n,int m>
	class inner_product_topology:public norm_topology<K,n,m>
	{
		using matrix_type = norm_topology<K, n,m>::matrix_type;
	public:
		virtual K inner_product(const matrix_type& p, const matrix_type& q) const = 0;
		virtual K dot_product(const matrix_type& p, const matrix_type& q) const
		{
			return inner_product(p.conj(), q);
		}
		virtual real_field norm(const matrix_type& p) const
		{
			return std::sqrt(static_cast<real_field>(inner_product(p, p)));
		}
	};

	template<typename K, int n,int m >
	class L2_vect_inner_product:public inner_product_topology<K,n,m>
	{
		using matrix_type = norm_topology<K, n,m>::matrix_type;
	public:
		K inner_product(const matrix_type& P, const matrix_type& Q) const
		{
			auto W = P.conj();
			K R;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					R += W[i][j] * Q[i][j];
			return R;
		}
	};
}