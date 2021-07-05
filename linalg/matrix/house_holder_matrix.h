#pragma once
#include "linalg/matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::linalg::special
{
	template<typename K, int n>
	class house_holder_matrix :public matrix<K, n,n>
	{
	public:
		house_holder_matrix(const finite_dimensional_vector_space<K, n>& v)
		{
			auto u=v / K(v.norm());
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					this->u[i][j] = K(i == j) - K(2)*u[i] * u[j].conj();

		}

		house_holder_matrix(const std::vector<K>& v):house_holder_matrix(finite_dimensional_vector_space<K,n>(v))
		{
	
		}
	};

}