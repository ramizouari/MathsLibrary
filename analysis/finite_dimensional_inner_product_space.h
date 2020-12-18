#pragma once
#include "linalg/finite_dimensional_vector_space.h"
#include "complex.h"
#include "normed_finite_dimensional_space.h"
#include "linalg/matrix.h"

template<int n>
class Lp_finite_dimensional_space <::complex, 2, n > :virtual public finite_dimensional_vector_space<::complex, n>,
	virtual public normed_space<::complex>
{
	public:
		using finite_dimensional_vector_space <::complex, n>::finite_dimensional_vector_space;
		Lp_finite_dimensional_space(const finite_dimensional_vector_space<::complex, n>& a) :
			finite_dimensional_vector_space<::complex, n>(a)
		{
		}

		Lp_finite_dimensional_space(finite_dimensional_vector_space<::complex, n>&& a) :
			finite_dimensional_vector_space<::complex, n>(std::move(a))
		{
		}
		real_field norm() const override
		{
			std::vector<real_field> S;
			return std::sqrt(
				std::reduce(this->u.begin(), this->u.end(), real_field(0), 
					[](const auto& a, const auto& b) {return a + std::pow(b.abs(),2); })
			);
		}

		::complex inner_product(const Lp_finite_dimensional_space& w) const
		{
			Lp_finite_dimensional_space s = w.conj();
			return std::inner_product(this->u.begin(), this->u.end(),s.u.begin(), ::complex(0));
		}
	};

template<typename F,int n>
using finte_dimensional_inner_product_space = Lp_finite_dimensional_space<F,2, n>;

template<int n>
using finite_dimensional_hilbert_space = Lp_finite_dimensional_space < ::complex, 2, n>;
template<int n>
using finite_dimensional_complex_space = finite_dimensional_hilbert_space<n>;



template<int n, int m>
finite_dimensional_hilbert_space<n> operator*(const matrix<::complex, n, m>& A,
	const finite_dimensional_hilbert_space<m>& u)
{
	finite_dimensional_hilbert_space<n> v;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			v.at(i) += A.at(i).at(j) * u.at(j);
	return v;
}

template<int n>
class Lp_finite_dimensional_space <real_field, 2, n > :virtual public finite_dimensional_vector_space<real_field, n>,
	virtual public normed_space<real_field>
{
public:
	using finite_dimensional_vector_space <real_field, n>::finite_dimensional_vector_space;
	Lp_finite_dimensional_space(const finite_dimensional_vector_space<real_field, n>& a) :
		finite_dimensional_vector_space<real_field, n>(a)
	{
	}

	Lp_finite_dimensional_space(finite_dimensional_vector_space<real_field, n>&& a) :
		finite_dimensional_vector_space<real_field, n>(std::move(a))
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