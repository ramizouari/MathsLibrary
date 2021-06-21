#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "prob/uniform_complex_generator.h"
#include "prob/uniform_real_generator.h"
#include "linalg/decomposer/QR_decomposition.h"

namespace math_rz::linalg::diagonalisation
{
	template<typename K, int n>
	class QR_algorithm
	{
		using eig_couple_vector = std::pair<coordinate_space<K, n>, square_matrix<K, n>>;
		using eig_couple = std::pair<K, coordinate_space<K, n>>;
		linalg::decomposer::QR_decomposition<K,n> QR;
		int steps = 50;
	public:
		finite_dimensional_vector_space<K,n> eigenvalues(const square_matrix<K, n>& _M) const 
		{
			
			auto M = _M;
			int p = steps;
			while (p--)
			{
				auto [Q, R] = QR.decompose(M);
				M = R * Q;
			}
			finite_dimensional_vector_space<K, n> E;
			for (int i = 0; i < n; i++)
				E[i] = M[i][i];
			return E;
		}
	
	};
}