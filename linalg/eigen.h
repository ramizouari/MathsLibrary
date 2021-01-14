#pragma once
#include "square_matrix.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "prob/uniform_real_generator.h"
namespace math_rz::linalg
{
	template<typename K,int n> requires vector_space_constraint::normed_vector_space<K>
	K largest_eig(const square_matrix<K, n>& M,const real_field &eps=1e-5)
	{
		if constexpr (n == 1)
			return 1;
		else {
			uniform_real_generator dist(-1, 1, 1000);
			analysis::Lp_finite_dimensional_space<K, 2, n> u, v, w;
			for (int i = 0; i < n; i++)
				u[i] = dist.generate();
			v = M * u;
			if (v.norm() == 0)
				return 0;
			v /= v.norm();
			while ((w = (v - u)).norm() > eps)
			{
				u = v;
				v = M * v;
				v /= v.norm();
			}
			int k = 0;
			for (int i = 0; i < n; i++)
				if (v[i].norm() > v[k].norm())
					k = i;
			if (v[k].norm() < eps)
				return 0;
			return (M * v)[k] / v[k];
		}
	}

	template<typename K,int n,int m> requires vector_space_constraint::normed_vector_space<K>
	real_field largest_sing(const matrix<K, n,m>& M, const real_field& eps = 1e-5)
	{
		return std::sqrt(static_cast<real_field>(largest_eig(square_matrix<K,m>(M.conj_transpose()*M))));
	}

	template<typename K, int n, int m>
	void gram_schmidt_inplace(matrix<K,n,m>&M )
	{
		std::vector<coordinate_space<K, m>> orthonormal_vectors;
		for (int i = 0; i < n; i++)
		{
			coordinate_space<K, m> w(M[i]);
			for (auto& s : orthonormal_vectors)
				w -= s.inner_product(w) * s;
			if(!w.is_zero())
				w /= w.norm();
			orthonormal_vectors.push_back(w);
			M[i] = orthonormal_vectors[i].get_vect();
		}
	}

	template<typename K, int n, int m>
	matrix<K,n,m> gram_schmidt(const matrix<K, n, m>& M)
	{
		matrix<K, n, m> W(M);
		gram_schmidt_inplace<K,n,m>(W);
		return W;
	}
}