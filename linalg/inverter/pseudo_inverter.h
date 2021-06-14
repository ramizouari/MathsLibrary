#pragma once
namespace math_rz::linalg::inverter
{
	template<typename K, int n, int m>
	class pseudo_inverter
	{
		using matrix_type = std::conditional_t<n == m, square_matrix<K, n>, matrix<K, n, m>>;
		using transpose_matrix_type = std::conditional_t<n == m, square_matrix<K, n>, matrix<K, m, n>>;
	public:
		virtual transpose_matrix_type pinv(const matrix_type& A) const = 0;
	};
}