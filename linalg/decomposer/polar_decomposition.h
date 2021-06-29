#pragma once
#include "decomposer.h"
#include "linalg/matrix/diagonal.h"
#include "linalg/diagonalisation/positive_semi_definite_diagonalisator.h"
#include "linalg/inverter/moore_penrose_pseudo_inverter.h"
#include "linalg/decomposer/SVD_decomposition.h"

namespace math_rz::linalg::decomposer
{
	template<typename K,int n,int m=n>
	struct UP
	{
		matrix<K, n, m> U;
		matrix<K, n, n> P;
	};
	template<typename K, int n, int m = n>
	class polar_decomposition :public decomposer<K,n,m,UP<K,n,m>>
	{
		SVD_decomposition<K, n, m> SVD;
	public:
		using UP = UP<K, n, m>;
		UP decompose(const matrix<K,n,m>&A) const
		{
			auto [U, D, V] = SVD.decompose(A);
			return UP{ U,D * V };
		}
	};
}