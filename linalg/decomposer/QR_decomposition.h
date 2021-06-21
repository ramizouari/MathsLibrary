#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "cholesky.h"
#include "decomposer.h"
namespace math_rz::linalg::decomposer
{
	template<typename K, int n, int m = n>
	struct QR
	{
		using matrix_type = std::conditional_t<n == m, square_matrix<K, n>, matrix<K, n, m>>;
		matrix_type Q;
		square_matrix<K, m> R;
	};
	template<typename K, int n,int m=n>
	class QR_decomposition:public decomposer<K,n,m,QR<K,n,m>>
	{
		using matrix_type = std::conditional_t<n == m, square_matrix<K, n>, matrix<K, n, m>>;
		using QR = QR<K, n, m>;
		cholesky<K,m> Ch;
	public:
		QR decompose(const matrix_type& A) const
		{
			auto R = Ch.orthogonalisation_cholesky(A.conj_transpose() * A).conj_transpose();
			auto Q = A*R.inv();
			return { Q,R };
		}
	};
}