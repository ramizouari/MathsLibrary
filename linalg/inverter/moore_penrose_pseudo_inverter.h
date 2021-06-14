#pragma once
#include "linalg/decomposer/QR_decomposition.h"
#include "pseudo_inverter.h"
namespace math_rz::linalg::inverter
{
	template<typename K, int n, int m>
	class moore_penrose_pseudo_inverter :public pseudo_inverter<K, n, m>
	{
		using matrix_type = std::conditional_t<n == m, square_matrix<K, n>, matrix<K, n, m>>;
		using transpose_matrix_type = std::conditional_t<n == m, square_matrix<K, n>, matrix<K, m, n>>;
	public:
		virtual transpose_matrix_type pinv(const matrix_type& A) const override
		{
			static math_rz::linalg::decomposer::QR_decomposition<K,n,m> QR;
			auto [Q, R] = QR.decompose(A);
			return R.inv() * Q.H();
		}
	};
}
