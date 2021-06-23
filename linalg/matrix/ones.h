#pragma once
#pragma once
#include "linalg/matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::linalg::special
{
	template<typename K, int n,int m>
	class ones :public matrix<K, n,m>
	{
	public:
		ones()
		{
			for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
				this->u[i][j] = K(1);
		}
	};

}