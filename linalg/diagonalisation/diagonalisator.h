#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "prob/uniform_complex_generator.h"
#include "prob/uniform_real_generator.h"

namespace math_rz::linalg::diagonalisation
{
	template<typename K,int n>
	struct eigenbasis
	{
		matrix<K, n, n> B;
		finite_dimensional_vector_space<K, n> D;
	};
	template<typename K,int n>
	class diagonalisator
	{
	public:
		using eigenbasis = eigenbasis<K, n>;
	};
}