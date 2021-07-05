#pragma once
#include "linalg/matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "complex.h"
#include <numbers>
namespace math_rz::linalg::special
{
	template<typename K, int n >
	class vandermonde_matrix :public matrix<K, n>
	{
	public:
		vandermonde_matrix(const finite_dimensional_vector_space<K,n> &_u)
		{
			for (int i = 0; i < n; i++)
			{
				K s = 1;
				for (int j = 0; j < n; j++) 
				{
					this->u[i][j] = s;
					s *= _u[i];
				}
			}

		}

		vandermonde_matrix(const std::vector<K>& _u) :vandermonde_matrix(finite_dimensional_vector_space<K, n>(_u))
		{

		}

	};

}