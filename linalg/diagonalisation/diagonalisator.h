#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "prob/uniform_complex_generator.h"
#include "prob/uniform_real_generator.h"
#include "linalg/decomposer/decomposer.h"

namespace math_rz::linalg::diagonalisation
{
	template<typename K,int n>
	struct eigenbasis
	{
		matrix<K, n, n> B;
		finite_dimensional_vector_space<K, n> D;
	};
	template<typename K,int n>
	class diagonalisator:public decomposer::decomposer<K,n,n,eigenbasis<K,n>>
	{
	public:
		using eigenbasis_type = eigenbasis<K, n>;
		virtual eigenbasis_type eigendecomposition(const matrix<K, n, n>& A) const = 0;
        eigenbasis_type decompose(const matrix<K, n, n>& A) const
		{
			return eigendecomposition(A);
		}

		virtual finite_dimensional_vector_space<K, n> eigenvalues(const matrix<K, n, n>& A) const
		{
			return eigendecomposition(A).D;
		}
	};
}