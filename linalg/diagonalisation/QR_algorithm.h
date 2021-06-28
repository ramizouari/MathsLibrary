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


namespace math_rz::linalg::diagonalisation
{
	template<typename K, int n>
	class QR_algorithm :public diagonalisator<K,n>
	{
		using eig_couple_vector = std::pair<coordinate_space<K, n>, matrix<K, n, n>>;
		using eig_couple = std::pair<K, coordinate_space<K, n>>;
		linalg::decomposer::QR_decomposition<K,n> QR;
		upper_triangular_diagonalisator<K, n> UTD;
		int steps = 50;
		real_field eps = 1e-5;

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
	public:
		struct eigenbasis
		{
			matrix<K, n, n> B;
			finite_dimensional_vector_space<K, n>D;
		};

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
			while (p--)
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
}