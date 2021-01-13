#pragma once
#include "linalg/finite_dimensional_vector_space.h"
#include "complex.h"
#include "normed_finite_dimensional_space.h"
#include "linalg/matrix.h"


namespace math_rz::analysis {
	template<int n>
	class Lp_finite_dimensional_space <complex, 2, n > :virtual public linalg::finite_dimensional_vector_space<complex, n>,
		virtual public normed_space<complex>
	{
	public:
		using linalg::finite_dimensional_vector_space <complex, n>::finite_dimensional_vector_space;
		Lp_finite_dimensional_space(const linalg::finite_dimensional_vector_space<complex, n>& a) :
			linalg::finite_dimensional_vector_space<complex, n>(a)
		{
		}

		Lp_finite_dimensional_space(linalg::finite_dimensional_vector_space<complex, n>&& a) :
			linalg::finite_dimensional_vector_space<complex, n>(std::move(a))
		{
		}
		real_field norm() const override
		{
			std::vector<real_field> S;
			return std::sqrt(
				std::reduce(this->u.begin(), this->u.end(), real_field(0),
					[](const auto& a, const auto& b) {return a + std::pow(b.abs(), 2); })
			);
		}

		complex inner_product(const Lp_finite_dimensional_space& w) const
		{
			Lp_finite_dimensional_space s = w.conj();
			return std::inner_product(this->u.begin(), this->u.end(), s.u.begin(), complex(0));
		}
	};

	template<typename K, int n>
	using finte_dimensional_inner_product_space = Lp_finite_dimensional_space<K, 2, n>;

	template<int n>
	using finite_dimensional_hilbert_space = Lp_finite_dimensional_space < complex, 2, n>;
	template<int n>
	using finite_dimensional_complex_space = finite_dimensional_hilbert_space<n>;



	

	template<int n>
	class Lp_finite_dimensional_space <real_field, 2, n > :virtual public linalg::finite_dimensional_vector_space<real_field, n>,
		virtual public normed_space<real_field>
	{
	public:
		using linalg::finite_dimensional_vector_space <real_field, n>::finite_dimensional_vector_space;
		Lp_finite_dimensional_space(const linalg::finite_dimensional_vector_space<real_field, n>& a) :
			linalg::finite_dimensional_vector_space<real_field, n>(a)
		{
		}

		Lp_finite_dimensional_space(linalg::finite_dimensional_vector_space<real_field, n>&& a) :
			linalg::finite_dimensional_vector_space<real_field, n>(std::move(a))
		{
		}
		real_field norm() const override
		{
			std::vector<real_field> S;
			return std::sqrt(
				inner_product(*this)
			);
		}

		real_field inner_product(const Lp_finite_dimensional_space& w) const
		{
			return std::inner_product(this->u.begin(), this->u.end(), w.u.begin(), real_field(0));
		}
	};

	template<int n>
	using euclidean_space = L2_finite_dimensional_space<real_field, n>;


	namespace finite_dimensional
	{
		template<int n>
		using hilbert_space = finite_dimensional_hilbert_space<n>;

		template<int n>
		using complex_space = finite_dimensional_complex_space<n>;

		template<int n>
		using euclidean_space = L2_finite_dimensional_space<real_field, n>;
	}
}