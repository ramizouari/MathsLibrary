#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "cholesky.h"
#include "decomposer.h"
#include <unordered_set>
namespace math_rz::linalg::decomposer
{
	template<typename K, int n, int m = n>
	struct QR
	{
		matrix<K,n,m> Q;
		matrix<K, m,m> R;
	};
	template<typename K, int n,int m=n>
	class QR_decomposition:public decomposer<K,n,m,QR<K,n,m>>
	{
		using matrix_type = matrix<K,n,m>;
		using QR_type = QR<K, n, m>;
	protected:
		cholesky<K,m> Ch;
		real_field eps = 1e-7;

		std::vector<finite_dimensional_vector_space<K, n>> gram_schmidt(const std::vector<finite_dimensional_vector_space<K, n>>& A)const
		{

			std::vector<finite_dimensional_vector_space<K, n>> B;
			std::unordered_set<int> J;
			finite_dimensional_vector_space<K, n> e;
			for (int i = 0; i < A.size(); i++)
			{
				finite_dimensional_vector_space<K, n> R, u;
				for (auto& s : B)
					R += A[i].inner_product(s) * s;
				u = A[i] - R;
				real_field N = u.norm();
				if (N > eps)
				{
					u /= N;
					J.insert(i);
					B.push_back(u);
				}
			}
			return B;
		}

		std::vector<finite_dimensional_vector_space<K, n>> gram_schmidt(const matrix<K, n, n>& A)const
		{
			std::vector<finite_dimensional_vector_space<K, n>> U(n);
			for (int i = 0; i < n; i++)
				U[i] = A.get_column(i);
			return gram_schmidt(U);
		}
	public:
        QR_type decompose(const matrix_type& A) const
		{
			auto R = Ch.iterative_cholesky(A.conj_transpose() * A).conj_transpose();
			std::vector<finite_dimensional_vector_space<K, n>> G;
			auto Q = A*R.inv();
			for (int i = 0; i < n; i++)
				G.push_back(Q.get_column(i));
			for (int i = 0; i < n; i++)
			{
				finite_dimensional_vector_space<K, n> e;
				e[i] = 1;
				G.push_back(e);
			}
			auto B= gram_schmidt(G);
			for (int i = 0; i < n; i++)
				for(int j = 0; j < n; j++)
					Q[i][j] = B[j][i];
			return { Q,R };
		}
	};
	template<typename K, int n, int m = n>
	class QR_decomposition_lite :public QR_decomposition<K,n,m>
	{
		using matrix_type = matrix<K, n, m>;
		using QR_type = QR<K, n, m>;

	public:
        QR_type decompose(const matrix_type& A) const
		{
			auto R = this->Ch.iterative_cholesky(A.conj_transpose() * A).conj_transpose();
			auto Q = A * R.inv();
			return { Q,R };
		}
	};
}