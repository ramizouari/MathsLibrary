#pragma once
#include "linalg/square_matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "prob/uniform_complex_generator.h"
#include "prob/uniform_real_generator.h"

namespace math_rz::linalg::diagonalisation
{
	template<typename K,int n>
	class diagonalisator
	{
		using eig_couple_vector = std::pair<coordinate_space<K, n>, matrix<K, n, n>>;
		using eig_couple = std::pair<K, coordinate_space<K, n>>;
	public:
		virtual eig_couple_vector diagonalise(const matrix<K, n, n>& M) const
		{

		}
	protected:
		virtual eig_couple eigen_couple(const matrix<K,n,n>&M) const
		{
			if constexpr (n == 1)
				return { 1,1 };
			else {
				coordinate_space<K, 2> u, v, w;
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
				for (int i = 0; i < n; i++)
					M[i][i] += delta;
				v = M * u;
				if (v.norm() == 0)
					return { delta,u };
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
					return { -delta,u };
				coordinate_space<K, n> s = M * v;
				return { (s[k] / v[k])-delta,s };
			}
		}
		virtual eig_couple_vector diagonalise_inplace(matrix<K,n,n> &M) const
		{
			eig_couple_vector C;
			K offset;
			for (int i = 0; i < n; i++)
			{
				auto [lambda,u]= eigen_couple(M);
				C.first[i] = lambda + offset;
				C.second[i] = u;

			}
		}
	};
}