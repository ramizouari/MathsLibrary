#pragma once
#include "normed_space.h"
#include <algorithm>
#include "linalg/finite_dimensional_vector_space.h"
#include <iterator>
#include "linalg/matrix.h"

namespace math_rz
{
	template<field_constraints::has_abs K, int p, int n >
	class Lp_finite_dimensional_space : virtual public normed_space<K>,
		virtual public finite_dimensional_vector_space<K, n>
	{
	public:
		using finite_dimensional_vector_space < K, n>::finite_dimensional_vector_space;



		Lp_finite_dimensional_space(const finite_dimensional_vector_space<K, n>& a) :
			finite_dimensional_vector_space<K, n>(a)
		{
		}

		Lp_finite_dimensional_space(finite_dimensional_vector_space<K, n>&& a) :
			finite_dimensional_vector_space<K, n>(std::move(a))
		{
		}
		real_field norm() const override
		{
			return
				std::pow(
					std::reduce(this->u.begin(), this->u.end(), real_field(0), [](const real_field& a, const K& b)
						{
							return a + std::pow(b.abs(), p); }),
					real_field(1. / p)
								);

		}
	};

	template<field_constraints::has_abs K, int n>
	class Linf_finite_dimensional_space : virtual public normed_space<K>,
		virtual public finite_dimensional_vector_space<K, n>
	{
	public:
		using finite_dimensional_vector_space < K, n>::finite_dimensional_vector_space;



		Linf_finite_dimensional_space(const finite_dimensional_vector_space<K, n>& a) :
			finite_dimensional_vector_space<K, n>(a)
		{
		}

		Linf_finite_dimensional_space(finite_dimensional_vector_space<K, n>&& a) :
			finite_dimensional_vector_space<K, n>(std::move(a))
		{
		}
		real_field norm() const override
		{
			return std::reduce(this->u.begin(), this->u.end(), real_field(0), [](const real_field& a, const K& b) {
				return std::max<real_field>(a, b.abs()); });

		}
	};

	template<int k>
	concept multi_dimensional = k > 1;

	template <typename K, int n, int m> requires multi_dimensional<n>
	Linf_finite_dimensional_space<K, n> operator*(const matrix<K, n, m>& A, const Linf_finite_dimensional_space<K, m>& u)
	{
		Linf_finite_dimensional_space<K, n> v;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				v.at(i) += A.at(i).at(j) * u.at(j);
		return v;
	}

	template <typename K, int p, int n, int m>requires multi_dimensional<n>
	Lp_finite_dimensional_space<K, p, n> operator*(const matrix<K, n, m>& A,
		const Lp_finite_dimensional_space<K, p, m>& u)
	{
		Lp_finite_dimensional_space<K, p, n> v;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				v.at(i) += A.at(i).at(j) * u.at(j);
		return v;
	}


	template <typename K, int p, int n>
	Lp_finite_dimensional_space<K, p, n> operator*(const K& r,
		const Lp_finite_dimensional_space<K, p, n>& u)
	{
		Lp_finite_dimensional_space<K, p, n> v(u);
		for (int i = 0; i < n; i++)
			v.at(i) *= r;
		return v;
	}

	template <typename K, int p, int n>
	Lp_finite_dimensional_space<K, p, n> operator*(const long double& r,
		const Lp_finite_dimensional_space<K, p, n>& u)
	{
		Lp_finite_dimensional_space<K, p, n> v(u);
		for (int i = 0; i < n; i++)
			v.at(i) *= r;
		return v;
	}

	template<typename K, int n >
	using L1_finite_dimensional_space = Lp_finite_dimensional_space<K, 1, n>;

	template<typename K, int n>
	using L2_finite_dimensional_space = Lp_finite_dimensional_space < K, 2, n>;


	namespace finite_dimensional
	{
		template<typename K, int n>
		using L1_space = L1_finite_dimensional_space<K, n>;

		template<typename K, int n>
		using L2_space = L2_finite_dimensional_space<K, n>;

		template<typename K, int p, int n >
		using Lp_space = Lp_finite_dimensional_space<K, p, n>;

	}
}