#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz::linalg::special
{
	template<typename K, int n>
	class givens_matrix :public matrix < K, n, n >
	{
	public:
		givens_matrix(int p,int q,K a,K b)
		{
			K r = std::sqrt(pow(a.norm(),2) + pow(b.norm(),2));
			a /= r;
			b /= r;
			for (int i = 0; i < n; i++)
				if (i != p && i != q)
					this->u[i][i] = 1;
			this->u[p][p] = a;
			this->u[q][p] = b;
			this->u[p][q] = -b;
			this->u[q][q] = a;


		}

		givens_matrix(int i, int j,real_field theta) :givens_matrix(i,j,std::cos(theta),std::sin(theta))
		{

		}
	};

}