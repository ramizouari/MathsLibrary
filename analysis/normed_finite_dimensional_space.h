#pragma once
#include "normed_space.h"
#include <algorithm>
#include "linalg/finite_dimensional_vector_space.h"
#include <iterator>
#include "linalg/matrix.h"

namespace math_rz
{
	template<field_constraints::has_abs F, int p, int n >
	class Lp_finite_dimensional_space : virtual public normed_space<F>,
		virtual public finite_dimensional_vector_space<F, n>
	{
	public:
		using finite_dimensional_vector_space < F, n>::finite_dimensional_vector_space;



		Lp_finite_dimensional_space(const finite_dimensional_vector_space<F, n>& a) :
			finite_dimensional_vector_space<F, n>(a)
		{
		}

		Lp_finite_dimensional_space(finite_dimensional_vector_space<F, n>&& a) :
			finite_dimensional_vector_space<F, n>(std::move(a))
		{
		}
		real_field norm() const override
		{
			return
				std::pow(
					std::reduce(this->u.begin(), this->u.end(), real_field(0), [](const real_field& a, const F& b)
						{
							return a + std::pow(b.abs(), p); }),
					real_field(1. / p)
								);

		}
	};

	template<field_constraints::has_abs F, int n>
	class Linf_finite_dimensional_space : virtual public normed_space<F>,
		virtual public finite_dimensional_vector_space<F, n>
	{
	public:
		using finite_dimensional_vector_space < F, n>::finite_dimensional_vector_space;



		Linf_finite_dimensional_space(const finite_dimensional_vector_space<F, n>& a) :
			finite_dimensional_vector_space<F, n>(a)
		{
		}

		Linf_finite_dimensional_space(finite_dimensional_vector_space<F, n>&& a) :
			finite_dimensional_vector_space<F, n>(std::move(a))
		{
		}
		real_field norm() const override
		{
			return std::reduce(this->u.begin(), this->u.end(), real_field(0), [](const real_field& a, const F& b) {
				return std::max<real_field>(a, b.abs()); });

		}
	};


	template <typename F, int n, int m>
	Linf_finite_dimensional_space<F, n> operator*(const matrix<F, n, m>& A, const Linf_finite_dimensional_space<F, m>& u)
	{
		Linf_finite_dimensional_space<F, n> v;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				v.at(i) += A.at(i).at(j) * u.at(j);
		return v;
	}

	template <typename F, int p, int n, int m>
	Lp_finite_dimensional_space<F, p, n> operator*(const matrix<F, n, m>& A,
		const Lp_finite_dimensional_space<F, p, m>& u)
	{
		Lp_finite_dimensional_space<F, p, n> v;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				v.at(i) += A.at(i).at(j) * u.at(j);
		return v;
	}

	template<typename F, int n >
	using L1_finite_dimensional_space = Lp_finite_dimensional_space<F, 1, n>;

	template<typename F, int n>
	using L2_finite_dimensional_space = Lp_finite_dimensional_space < F, 2, n>;


	namespace finite_dimensional
	{
		template<typename F, int n>
		using L1_space = L1_finite_dimensional_space<F, n>;

		template<typename F, int n>
		using L2_space = L2_finite_dimensional_space<F, n>;

		template<typename F, int p, int n >
		using Lp_space = Lp_finite_dimensional_space<F, p, n>;

	}
}