#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::linalg::special
{
	template<typename K, int n>
	class circulant_matrix :public square_matrix<K, n>
	{
	public:
		circulant_matrix(const finite_dimensional_vector_space<K, n>& v)
		{
			auto u = v / K(v.norm());
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					this->u[i][j] = v[(j+n-i)%n];

		}

		circulant_matrix(const std::vector<K>& v) :circulant_matrix(finite_dimensional_vector_space<K, n>(v))
		{

		}
	};

}