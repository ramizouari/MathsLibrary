#pragma once
#include "diagonalisator.h"
#include "linalg/inverter/moore_penrose_pseudo_inverter.h"
#include <unordered_map>

namespace math_rz::linalg::diagonalisation
{
	template<typename K, int n>
	class upper_triangular_diagonalisator :public diagonalisator<K, n>
	{
		real_field eps;
		std::vector<finite_dimensional_vector_space<K, n>> gram_schmidt(const std::vector<finite_dimensional_vector_space<K, n>>& A)const
		{

			std::vector<finite_dimensional_vector_space<K, n>> B;
			std::unordered_set<int> J;
			finite_dimensional_vector_space<K, n> e;
			for (int i = 0; i < n; i++)
			{
				finite_dimensional_vector_space<K, n> R, u;
				for (auto& s : B)
					R += s.inner_product(A[i]) * s;
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
		inverter::moore_penrose_pseudo_inverter<K, n, n> PI;

	public:
		using eigenbasis=diagonalisator<K, n>::eigenbasis;
		
		upper_triangular_diagonalisator(real_field _eps = 1e-5) :eps(_eps) {}
		eigenbasis diagonalise(const matrix<K, n, n>& A) const
		{
			auto [B, N] = gram_schmidt(A);
			eigenbasis EB;
			EB.B = std::move(B);
			linalg::structure::vector::L2_induced_vect_inner_product<K, n> P(A);
			for (int i = 0; i < n; i++)
			{
				finite_dimensional_vector_space<K, n> u(EB.B.get_column(i));
				EB.D[i] = P.inner_product(u, u);
			}
			return EB;
		}

		finite_dimensional_vector_space<K, n> eigenvalues(const matrix<K, n, n>& _M) const
		{
			finite_dimensional_vector_space<K, n> u;
			for (int i = 0; i < n; i++)
				u[i] = _M[i][i];
			return u;
		}

		eigenbasis eigendecomposition(const matrix<K, n, n>& _M) const
		{
			std::vector<finite_dimensional_vector_space<K, n>> B;
			matrix<K, n, n> E;
			matrix<K, n, n> M(_M);
			std::vector <std::pair < K, int> > R;
			std::unordered_map<int, int> mapper;
			auto L = eigenvalues(M);
			/* Counting eigenvalues each with its correspending multiplicity*/
			for (int i = 0; i < n; i++)
			{
				bool different = true;
				for (int j = 0; j < i; j++)
					if (L[j].metric(L[i]) <= eps)
					{
						mapper[i] = mapper[j];
						R[mapper[i]].second++;
						different = false;
						break;
					}
				if (different)
				{
					mapper.insert({ i, R.size() });
					R.push_back({ L[i],1 });
				}

			}
			int s = 0;
			/*
			* Finding each eigenvalue
			*/
			for (auto I : R)
			{
				auto [lambda, m] = I;
				/*
				* Creating The matrix M-lambda*I
				*/
				for (int r = 0; r < n; r++)
					M[r][r] -= lambda;
				/*
				* Raising it to the nth power to detect potential generalized eigenvalues
				*/
				auto N = pow(M, 1);

				/*
				* Finding an orthonormal basis of the correspending eigenspace
				*/
				std::cout << N;
				auto C = matrix<K, n, n>(1) - PI.pinv(N) * N;
				std::cout << C;
				auto B1 = gram_schmidt(C);
				for (int a = 0; a < m; a++)
					B.push_back(B1[a]);
				for (int r = 0; r < n; r++)
					M[r][r] += lambda;
				s += m;
			}
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					E[i][j] = B[j][i];
			eigenbasis EB = { E,L };
			return EB;
		}
	};
}