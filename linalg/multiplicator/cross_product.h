#pragma once
#include "linalg/finite_dimensional_vector_space.h"
#include "analysis/finite_dimensional_inner_product_space.h"
namespace math_rz::analysis
{
	template<typename K>
	class cross_product
	{
	public:
		L2_finite_dimensional_space<K, 3> 
			multiply(const L2_finite_dimensional_space<K, 3>& u, 
				const L2_finite_dimensional_space<K, 3>& v) const
		{
			return L2_finite_dimensional_space<K, 3>({
				u[1] * v[2] - u[2] * v[1],
				u[2] * v[0] - u[0] * v[2],
				u[0] * v[1] - u[1] * v[0]
				});
		}
	};
}