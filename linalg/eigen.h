#pragma once
#include "square_matrix.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "prob/uniform_real_generator.h"
#include "prob/uniform_complex_generator.h"
namespace math_rz::linalg
{
	template<typename K, int n> requires vector_space_constraint::normed_vector_space<K>
		K largest_eig(const square_matrix<K, n>& M, const real_field& eps = 1e-5)
		{
			if constexpr (n == 1)
				return 1;
			else {
				uniform_real_generator dist(-1, 1, 1000);
				analysis::Lp_finite_dimensional_space<K, 2, n> u, v, w;
				real_field a, b;
				for (int i = 0; i < n; i++)
					u[i] = dist.generate();
				v = M * u;
				if (v.norm() == 0)
					return 0;
				v /= v.norm();
				do
				{
					u = v;
					v = M * v;
					a = b;
					b = v.norm() / u.norm();
					v /= v.norm();
				} while (std::abs(b - a) > eps);
				int k = 0;
				for (int i = 0; i < n; i++)
					if (v[i].norm() > v[k].norm())
						k = i;
				if (v[k].norm() < eps)
					return 0;
				return (M * v)[k] / v[k];
			}
		}

	template<typename K,int n> requires vector_space_constraint::normed_vector_space<K>
	K largest_pos_eig(const square_matrix<K, n>& M,const real_field &eps=1e-5)
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

	template<typename K, int n> requires vector_space_constraint::normed_vector_space<K>
		std::pair<K, coordinate_space<K, n>> 
			largest_eig_couple(const square_matrix<K, n>& M, const real_field& eps = 1e-5,
				const std::vector < std::pair<K, coordinate_space<K, n>>> P = {})
		{
			if constexpr (n == 1)
				return { 1,1 };
			else {
				analysis::Lp_finite_dimensional_space<K, 2, n> u, v, w;
				real_field a, b;
				K delta;
				if constexpr (std::is_same_v<K, complex>)
				{
					uniform_complex_generator dist(-1, 1, -1, 1, 1000);
					for (int i = 0; i < n; i++)
						u[i] = dist.generate();
					delta = dist.generate();
				}
				else
				{
					uniform_real_generator dist(-1, 1, 1000);
					for (int i = 0; i < n; i++)
						u[i] = dist.generate();
					delta = dist.generate();
				}
				v = M * u;
				if (v.norm() == 0)
					return { 0,u };
				v /= v.norm();
				do
				{
					u = v;
					v = M * v;
					for (auto [t, u] : P)
						v -=  u.inner_product(v)*u;
					a = b;
					b = v.norm() / u.norm();
					v /= v.norm();
				} while (std::abs(b - a) > eps);
				int k = 0;
				for (int i = 0; i < n; i++)
					if (v[i].norm() > v[k].norm())
						k = i;
				if (v[k].norm() < eps)
					return { 0,u };
				coordinate_space<K, n> s = M * v;
				for (auto [t, u] : P)
					s -=  u.inner_product(s) * u;
				return { (s[k] / v[k]),s };
			}
		}

	template<typename K, int n> requires vector_space_constraint::normed_vector_space<K>
	std::pair<K,coordinate_space<K,n>> 
		largest_pos_eig_couple(const square_matrix<K, n>& M, const real_field& eps = 1e-5)
		{
			if constexpr (n == 1)
				return { 1,1 };
			else {
				analysis::Lp_finite_dimensional_space<K, 2, n> u, v, w;
				if constexpr (std::is_same_v<K,complex>)
				{ 
					uniform_complex_generator dist(-1, 1,-1,1, 1000);
					for (int i = 0; i < n; i++)
						u[i] = dist.generate();
				}
				else
				{
				uniform_real_generator dist(-1, 1, 1000);
					for (int i = 0; i < n; i++)
						u[i] = dist.generate();
				}
				v = M * u;
				if (v.norm() == 0)
					return { 0,u };
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
					return { 0,u };
				coordinate_space<K,n> s = M * v;
				return { (s[k] / v[k]),s};
			}
		}

	template<typename K,int n,int m> requires vector_space_constraint::normed_vector_space<K>
	real_field largest_sing(const matrix<K, n,m>& M, const real_field& eps = 1e-5)
	{
		return std::sqrt(static_cast<real_field>(largest_pos_eig(square_matrix<K,m>(M.conj_transpose()*M))));
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

	template<typename K,int n>
	square_matrix<K, n> proj(const coordinate_space<K, n>& w)
	{
		square_matrix<K, n> M;
		coordinate_space<K, n> s;
		for (int i = 0; i < n; i++)
		{
			s = w[i] * w;
			M[i] = s.get_vect();
		}
		return M.T();
	}

	template<typename K, int n> 
	requires vector_space_constraint::normed_vector_space<K>
	std::pair<coordinate_space<K,n>, square_matrix<K, n>>
		eigdecomposition_inplace(square_matrix<K, n>& M, const real_field& eps = 1e-5)
	{
		std::vector<std::vector<K>> U(n,std::vector<K>(n));
		std::vector<K> D(n);
		std::vector<std::pair<K, coordinate_space<K, n>>> P;
		for (int i = 0; i < n; i++)
		{
			auto [t, u] = largest_eig_couple(M, eps,P);
			auto N = u.norm();
			if (N == 0)
				break;
			u /= N;
			U[i] = u.get_vect();
			D[i] = t;
			auto S=t*proj(u);
			P.push_back({ D[i],U[i] });
			//M -= S;
		}
		return { coordinate_space<K,n>(D),square_matrix<K,n>(U).T() };
	}

	template<typename K, int n> requires vector_space_constraint::normed_vector_space<K>
		std::pair<coordinate_space<K, n>, square_matrix<K, n>>
			eigdecomposition(const square_matrix<K, n>& M, const real_field& eps = 1e-5)
		{
			square_matrix<K, n> S(M);
			return eigdecomposition_inplace<K, n>(S, eps);
		}
}