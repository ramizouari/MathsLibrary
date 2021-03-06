#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "cholesky.h"
#include "decomposer.h"
#include "QR_decomposition.h"
#include "LQ_decomposition.h"
#include "linalg/diagonalisation/QR_algorithm.h"
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
		using SVD_type = SVD<K, n, m>;
		QR_decomposition<K, n, m> QR;
		LQ_decomposition<K, n, m> LQ;
		integer steps = 100;
		cholesky<K, m> Ch;
		diagonalisation::QR_algorithm_hermitian < K, std::min(n, m)> QR_diag;
	public:
        SVD_type decompose(const matrix_type& A) const
		{
            SVD_type UDV;
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

        SVD_type svd_decomposition(const matrix_type& A) const
		{
            return decompose(A);
		}

		finite_dimensional_vector_space<K, std::min(m,n)> singular_values(const matrix_type&A) const
		{
			if constexpr (n <= m)
				return QR_diag.eigenvalues(A * A.H());
			else return QR_diag.eigenvalues(A.H() * A);
		}
	};
}