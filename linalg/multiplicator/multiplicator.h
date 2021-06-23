#pragma once
#include "linalg/matrix.h"
#include <type_traits>
#include <bit>
#include <thread>
#include <future>
#include <atomic>

namespace math_rz::linalg
{
	template<typename K>
	class multiplicator
	{
		template<int n,int m>
		using matrix_type = std::conditional_t<n == m, matrix<K, n,n>, matrix<K, n, m>>;
	public:
		template<int n,int m,int p>
		matrix_type<n, p> multiply(const matrix_type<n, m>& A, const matrix_type<m, p>& B) const
		{
			matrix_type<n, p> C;
//#pragma omp parallel for
			for (int i = 0; i < n; i++)
				for (int k = 0; k < p; k++)
					for (int j = 0; j < m; j++)
						C[i][j] += A[i][k] * B[k][j];
			return C;
		}
	};

	constexpr unsigned int bit_ceil(unsigned int s)
	{
		int r = 1;
		while (r < s) r *= 2;
		return r;
	}


	template<typename K>
	class strassen_multiplicator : public multiplicator<K>
	{
		template<int n, int m>
		using matrix_type = std::conditional_t<n == m, matrix<K, n, n>, matrix<K, n, m>>;
	protected:
		static constexpr int size_limit = 64;
	public:
		template<int n, int m, int p>
		matrix_type<n, p> multiply(const matrix_type<n, m>& A, const matrix_type<m, p>& B) const
		{
			constexpr auto r = std::max({ n,m,p });
			constexpr int S = r+r%2;
			if (r <= size_limit)
				return multiplicator<K>::multiply<n, m, p>(A, B);
			matrix_type<S / 2, S / 2> A1, A2, A3, A4, B1, B2, B3, B4;
			for (int i = 0; i < std::min(S / 2, n); i++)
				for (int j = 0; j < std::min(S / 2, m); j++)
					A1[i][j] = A[i][j];

			for (int i = 0; i < std::min(S / 2, n); i++)
				for (int j = std::min(S / 2, m); j < m; j++)
					A2[i][j - S / 2] = A[i][j];

			for (int i = std::min(S / 2, n); i < n; i++)
				for (int j = 0; j < std::min(S / 2, m); j++)
					A3[i - S / 2][j] = A[i][j];

			for (int i = std::min(S / 2, n); i < n; i++)
				for (int j = std::min(S / 2, m); j < m; j++)
					A4[i - S / 2][j - S / 2] = A[i][j];

			for (int i = 0; i < std::min(S / 2, m); i++)
				for (int j = 0; j < std::min(S / 2, p); j++)
					B1[i][j] = B[i][j];

			for (int i = 0; i < std::min(S / 2, m); i++)
				for (int j = std::min(S / 2, p); j < p; j++)
					B2[i][j - S / 2] = B[i][j];

			for (int i = std::min(S / 2, m); i < m; i++)
				for (int j = 0; j < std::min(S / 2, p); j++)
					B3[i - S / 2][j] = B[i][j];

			for (int i = std::min(S / 2, m); i < m; i++)
				for (int j = std::min(S / 2, p); j < p; j++)
					B4[i - S / 2][j - S / 2] = B[i][j];
			matrix<K,S / 2,S/2> T1 = this->multiply<S / 2, S / 2,S/2>(A1 + A4, B1 + B4),
				T2 = this->multiply<S / 2, S / 2,S/2>(A3 + A4, B1),
				T3 = this->multiply<S / 2, S / 2,S/2>(A1, B2- B4),
				T4 = this->multiply<S / 2, S / 2,S/2>(A4, B3 - B1),
				T5 = this->multiply<S / 2, S / 2,S/2>(A1 + A2, B4),
				T6 = this->multiply<S / 2, S / 2,S/2>(A3 - A1, B1 + B2),
				T7 = this->multiply<S / 2, S / 2,S/2>(A2 - A4, B3 + B4);
			matrix_type<n, p> C;
			for (int i = 0; i < std::min(n, S / 2); i++)
				for (int j = 0; j < std::min(p, S / 2); j++)
					C[i][j] = T1[i][j] + T4[i][j] - T5[i][j] + T7[i][j];

			for (int i = 0; i < std::min(n, S / 2); i++)
				for (int j = std::min(p, S / 2); j < p; j++)
					C[i][j] = T3[i][j - S / 2] + T5[i][j - S / 2];

			for (int i = std::min(n, S / 2); i < n; i++)
				for (int j = 0; j < std::min(p, S / 2); j++)
					C[i][j] = T2[i - S / 2][j] + T4[i - S / 2][j];

			for (int i = std::min(n, S / 2); i < n; i++)
				for (int j = std::min(p, S / 2); j < p; j++)
					C[i][j] = T1[i - S / 2][j - S / 2] - T2[i - S / 2][j - S / 2] + T3[i - S / 2][j - S / 2] + T6[i - S / 2][j - S / 2];
			return C;
		}
	};

	template<typename K>
	class parallel_strassen_multiplicator : public strassen_multiplicator<K>
	{
		template<int n, int m>
		using matrix_type = std::conditional_t<n == m, matrix<K,n, n>, matrix<K, n, m>>;
		mutable std::atomic<int> threads = 0;
		int coeff;
	public:
		parallel_strassen_multiplicator(int _coeff = 1) :coeff(_coeff) {}
		template<int n, int m, int p>
		matrix_type<n, p> multiply(const matrix_type<n, m>& A, const matrix_type<m, p>& B) const
		{
			threads++;
			constexpr auto r = std::max({ n,m,p });
			if (threads > coeff * std::thread::hardware_concurrency() || r <= strassen_multiplicator<K>::size_limit)
			{
				threads--;
				return strassen_multiplicator<K>::multiply<n, m, p>(A, B);
			}
			constexpr int S = r+r%2;
			matrix_type<S / 2, S / 2> A1, A2, A3, A4, B1, B2, B3, B4;
			for (int i = 0; i < std::min(S / 2, n); i++)
				for (int j = 0; j < std::min(S / 2, m); j++)
					A1[i][j] = A[i][j];

			for (int i = 0; i < std::min(S / 2, n); i++)
				for (int j = std::min(S / 2, m); j < m; j++)
					A2[i][j - S / 2] = A[i][j];

			for (int i = std::min(S / 2, n); i < n; i++)
				for (int j = 0; j < std::min(S / 2, m); j++)
					A3[i - S / 2][j] = A[i][j];

			for (int i = std::min(S / 2, n); i < n; i++)
				for (int j = std::min(S / 2, m); j < m; j++)
					A4[i - S / 2][j - S / 2] = A[i][j];

			for (int i = 0; i < std::min(S / 2, m); i++)
				for (int j = 0; j < std::min(S / 2, p); j++)
					B1[i][j] = B[i][j];

			for (int i = 0; i < std::min(S / 2, m); i++)
				for (int j = std::min(S / 2, p); j < p; j++)
					B2[i][j - S / 2] = B[i][j];

			for (int i = std::min(S / 2, m); i < m; i++)
				for (int j = 0; j < std::min(S / 2, p); j++)
					B3[i - S / 2][j] = B[i][j];

			for (int i = std::min(S / 2, m); i < m; i++)
				for (int j = std::min(S / 2, p); j < p; j++)
					B4[i - S / 2][j - S / 2] = B[i][j];
			auto F1 = std::async(&parallel_strassen_multiplicator::multiply<S / 2, S / 2, S / 2>, this, A1 + A4, B1 + B4);
			auto F2 = std::async(&parallel_strassen_multiplicator::multiply<S / 2, S / 2, S / 2>, this, A3 + A4, B1);
			auto F3 = std::async(&parallel_strassen_multiplicator::multiply<S / 2, S / 2, S / 2>, this, A1, (B2 - B4));
			auto F4 = std::async(&parallel_strassen_multiplicator::multiply<S / 2, S / 2, S / 2>, this, (A4), B3 - B1);
			auto F5 = std::async(&parallel_strassen_multiplicator::multiply<S / 2, S / 2, S / 2>, this, A1 + A2, (B4));
			auto F6 = std::async(&parallel_strassen_multiplicator::multiply<S / 2, S / 2, S / 2>, this, A3 - A1, B1 + B2);
			auto F7 = std::async(&parallel_strassen_multiplicator::multiply<S / 2, S / 2, S / 2>, this, A2 - A4, B3 + B4);
			auto T1 = (F1.get());
			auto T2 = (F2.get());
			auto T3 = (F3.get());
			auto T4 = (F4.get());
			auto T5 = (F5.get());
			auto T6 = (F6.get());
			auto T7 = (F7.get());
			matrix_type<n, p> C;
			for (int i = 0; i < std::min(n, S / 2); i++)
				for (int j = 0; j < std::min(p, S / 2); j++)
					C[i][j] = T1[i][j] + T4[i][j] - T5[i][j] + T7[i][j];

			for (int i = 0; i < std::min(n, S / 2); i++)
				for (int j = std::min(p, S / 2); j < p; j++)
					C[i][j] = T3[i][j - S / 2] + T5[i][j - S / 2];

			for (int i = std::min(n, S / 2); i < n; i++)
				for (int j = 0; j < std::min(p, S / 2); j++)
					C[i][j] = T2[i - S / 2][j] + T4[i - S / 2][j];

			for (int i = std::min(n, S / 2); i < n; i++)
				for (int j = std::min(p, S / 2); j < p; j++)
					C[i][j] = T1[i - S / 2][j - S / 2] - T2[i - S / 2][j - S / 2] + T3[i - S / 2][j - S / 2] + T6[i - S / 2][j - S / 2];

			threads--;
			return C;
		}
	};
}