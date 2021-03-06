#pragma once
#include "linalg/matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::linalg::special
{
	template<typename K, int n>
	class centering_matrix :public matrix<K, n, n>
	{
	public:
		centering_matrix()
		{
			for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
				this->u[i][j] = K((i==j)-1./n);
		}
	};

}