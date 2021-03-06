#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/structure/vector/inner_product.h"
#include "decomposer.h"

namespace math_rz::linalg::decomposer
{
	template<typename K, int n>
	class cholesky : public decomposer < K, n, n, matrix < K, n,n >>
	{
		struct orthogonal_basis
		{
			matrix<K, n, n> B;
			std::vector<bool> N = std::vector<bool>(n, false);
		};
		real_field eps = 1e-8;
		
		orthogonal_basis gram_schmidt(const matrix<K,n,n>&A)const
		{

			std::vector<finite_dimensional_vector_space<K, n>> B(n);
			orthogonal_basis S;
			finite_dimensional_vector_space<K, n> e;
			linalg::structure::vector::L2_induced_vect_inner_product<K,n> P(A);
			for (int i = 0; i < n; i++)
			{
				e[i] =1;
				finite_dimensional_vector_space<K, n> R;
				for (int j = 0; j < i; j++)
					R += P.inner_product(B[j],e) * B[j];
				B[i] = e - R;
				real_field N = P.norm(B[i]);
				if (N > eps)
					B[i] /=  N;
				else S.N[i] = true;
				e[i] = 0;
			}
			for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
				S.B[i][j] = B[i][j];
			return S;
		}
	public:
		matrix<K, n ,n> iterative_cholesky(const matrix<K, n ,n>& A) const
		{
			matrix<K, n ,n> B;
			for (int j = 0; j < n; j++)
			{
				K S=A[j][j];
				for (int k = 0; k < j; k++)
					S -= B[j][k] * B[j][k].conj();
				B[j][j] = sqrt(S.abs());
				if (B[j][j].abs() < eps) for (int i = j + 1; i < n; i++)
					B[i][j] = 0;
				else for (int i = j+1; i < n; i++)
				{
					K S = A[i][j];
					for (int k = 0; k < j; k++)
						S -= B[i][k] * B[j][k].conj();
					B[i][j] = S / B[j][j];
				}
			}
			return B;
		}
		matrix<K, n ,n> decompose(const matrix<K, n ,n>& A) const override
		{
			return iterative_cholesky(A);
		}

		matrix<K, n ,n> orthogonalisation_cholesky(const matrix<K, n ,n>& A) const
		{
			orthogonal_basis S = gram_schmidt(A);
					matrix<K, n,n> D;
			for (int i = 0; i < n; i++) if (!S.N[i]) D[i][i] = 1;
			return  S.B.inv()*D;
		}
	};
}