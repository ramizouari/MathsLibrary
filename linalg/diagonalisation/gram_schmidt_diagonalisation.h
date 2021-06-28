#pragma once
#include "diagonalisator.h"
#include "positive_semi_definite_diagonalisator.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/matrix.h"
#include <vector>
namespace math_rz::linalg::diagonalisation
{
	template<typename K, int n>
	class gram_schmidt_diagonalisation :public positive_semi_definite_diagonalisator<K, n>
	{
		struct orthogonal_basis
		{
			matrix<K, n, n> B;
			std::vector<bool> N = std::vector<bool>(n, false);
		};
		real_field eps = 1e-8;

		orthogonal_basis gram_schmidt(const matrix<K, n, n>& A)const
		{

			std::vector<finite_dimensional_vector_space<K, n>> B(n);
			orthogonal_basis S;
			finite_dimensional_vector_space<K, n> e;
			linalg::structure::vector::L2_induced_vect_inner_product<K, n> P(A);
			for (int i = 0; i < n; i++)
			{
				e[i] = 1;
				finite_dimensional_vector_space<K, n> R;
				for (int j = 0; j < i; j++) if(!S.N[j])
					R += P.inner_product(e, B[j]) * B[j]/P.inner_product(B[j],B[j]);
				B[i] = e - R;
				real_field N = P.norm(B[i]);
				if (N <= eps)
					S.N[i] = true;
				e[i] = 0;
			}
			for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
				S.B[j][i] = B[i][j];
			return S;
		}

	public:
		struct eigenbase
		{
			matrix<K, n, n> B;
			finite_dimensional_vector_space<K, n> L;
		};
		gram_schmidt_diagonalisation(real_field _eps=1e-5):eps(_eps){}
		eigenbase diagonalise(const matrix<K,n,n> &A) const
		{
			auto [B, N] = gram_schmidt(A);
			eigenbase EB;
			EB.B=std::move(B);
			linalg::structure::vector::L2_induced_vect_inner_product<K, n> P(A);
			for (int i = 0; i < n; i++)
			{
				finite_dimensional_vector_space<K,n> u(EB.B.get_column(i));
				EB.L[i] = P.inner_product(u,u);
			}
			return EB;
		}
	};
}