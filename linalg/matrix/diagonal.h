#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::linalg::special
{
	template<typename K, int n>
	class diagonal :public matrix < K, n, n >
	{
	public:
		diagonal(const finite_dimensional_vector_space<K, n>& v)
		{
			for (int i = 0; i < n; i++)
				this->u[i][i] = v[i];
		}

		diagonal(const std::vector<K>& v) :diagonal(finite_dimensional_vector_space<K, n>(v))
		{

		}
	};

}