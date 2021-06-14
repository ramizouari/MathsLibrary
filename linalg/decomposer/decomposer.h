#pragma once
namespace math_rz::linalg::decomposer
{
	template<typename K, int n, int m, typename R>
	class decomposer
	{
		virtual R decompose(const matrix<K, n, m>& A) const = 0;
	};

	template<typename K, int n, typename R>
	class decomposer<K, n, n, R>
	{
		virtual R decompose(const square_matrix<K, n>& A) const = 0;
	};
}