#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "prob/uniform_complex_generator.h"
#include "prob/uniform_real_generator.h"
#include "linalg/decomposer/QR_decomposition.h"
#include "linalg/inverter/moore_penrose_pseudo_inverter.h"
#include <unordered_set>
#include <unordered_map>
#include "upper_triangular_diagonalisator.h"
#include "self_adjoint_diagonalisator.h"
#include "positive_semi_definite_diagonalisator.h"


namespace math_rz::linalg::diagonalisation
{
	template<typename K, int n>
	class QR_algorithm :virtual public diagonalisator<K,n>
	{
	protected:
		using eig_couple_vector = std::pair<coordinate_space<K, n>, matrix<K, n, n>>;
		using eig_couple = std::pair<K, coordinate_space<K, n>>;
		linalg::decomposer::QR_decomposition<K,n> QR;
		upper_triangular_diagonalisator<K, n> UTD;
		int steps = 200;
		real_field eps = 1e-7;

		std::vector<finite_dimensional_vector_space<K, n>> gram_schmidt(const std::vector<finite_dimensional_vector_space<K, n>>& A)const
		{

			std::vector<finite_dimensional_vector_space<K, n>> B;
			std::unordered_set<int> J;
			finite_dimensional_vector_space<K, n> e;
			for (int i = 0; i < n; i++)
			{
				finite_dimensional_vector_space<K, n> R, u;
				for (auto& s : B)
					R += A[i].inner_product(s)*s;
				u = A[i] - R;
				real_field N =u.norm();
				if(N > eps)
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

		real_field lower_triangular_norm(const matrix<K, n, n>& A) const
		{
			real_field S=0;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < i; j++)
					S = std::max(S, A[i][j].norm());
			return S;
		}

	public:
		using eigenbasis = diagonalisator<K, n>::eigenbasis;

		finite_dimensional_vector_space<K,n> eigenvalues(const matrix<K, n ,n>& _M) const 
		{
			
			auto M = _M;
			int p = steps;
			while (p--)
			{
				auto [Q, R] = QR.decompose(M);
				M = R * Q;
			}
			finite_dimensional_vector_space<K, n> E;
			for (int i = 0; i < n; i++)
				E[i] = M[i][i];
			return E;
		}

		eigenbasis eigendecomposition(const matrix<K, n, n>& _M) const
		{
			eigenbasis EB;
			EB.B = matrix<K, n, n>(1);
			auto M = _M;
			int p = steps;
			while (lower_triangular_norm(M) > eps &&p--)
			{
				auto [Q, R] = QR.decompose(M);
				M = R * Q;
				EB.B = EB.B*Q;
			}
			auto [V, D] = UTD.eigendecomposition(M);
			EB.B = EB.B*V;
			EB.D = D;
			return EB;
		}
	
	};

	template<typename K, int n>
	class QR_algorithm_hermitian :virtual public QR_algorithm<K,n>,virtual public self_adjoint_diagonalisator<K,n>
	{
	public:
		using eigenbasis = QR_algorithm<K, n>::eigenbasis;
		eigenbasis eigendecomposition(const matrix<K, n, n>& _M) const
		{
			eigenbasis EB;
			EB.B = matrix<K, n, n>(1);
			auto M = _M;
			int p = this->steps;
			while (this->lower_triangular_norm(M)>this->eps && p--)
			{
				auto [Q, R] = this->QR.decompose(M);
				M = R * Q;
				EB.B = EB.B * Q;
			}
			for (int i = 0; i < n; i++)
				EB.D[i] = M[i][i];
			return EB;
		}
	};

	template<typename K, int n>
	class QR_algorithm_positive_semi_definite :virtual public QR_algorithm_hermitian<K, n>, virtual public positive_semi_definite_diagonalisator<K, n>
	{
	};
}