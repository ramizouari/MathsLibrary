#pragma once
#include "linalg/matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "complex.h"
#include <numbers>
namespace math_rz::linalg::special
{
	template<int n>
	class DFT_matrix :public matrix<complex, n>
	{
		inline static complex w = std::exp(complex(0,-2* std::numbers::pi/n));
	public:
		DFT_matrix()
		{
			for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
			{
				this->u[i][j] = pow(static_cast<std::complex<long double>>(w), (i * j) % n);
				this->u[i][j] /= std::sqrt(n);
			}
			
		}
	};

}