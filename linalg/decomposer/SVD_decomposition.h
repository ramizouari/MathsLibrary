#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "cholesky.h"
#include "decomposer.h"
#include "QR_decomposition.h"
#include "LQ_decomposition.h"
namespace math_rz::linalg::decomposer
{
	template<typename K, int n, int m = n>
	struct SVD
	{
		matrix<K, n, n> U;
		matrix<K, n, m> D;
		matrix<K, m, m> V;
	};
	template<typename K, int n, int m = n>
	class SVD_decomposition :public decomposer<K, n, m, SVD<K, n, m>>
	{
		using matrix_type = matrix<K, n, m>;
		using SVD = SVD<K, n, m>;
		QR_decomposition<K, n, m> QR;
		LQ_decomposition<K, n, m> LQ;
		integer steps = 100;
		cholesky<K, m> Ch;
	public:
		SVD decompose(const matrix_type& A) const
		{
			SVD UDV;
			UDV.U = matrix<K, n, n>(1);
			UDV.V = matrix<K, n, n>(1);
			auto M = A;
			int p = steps;
			while (p--)
			{
				auto [Q, R] = QR.decompose(M);
				auto [L, P] = LQ.decompose(R);
				M = L;
				UDV.U = UDV.U * Q;
				UDV.V = P * UDV.V;
			}
			UDV.D = M;
			return UDV;
		}
	};
}