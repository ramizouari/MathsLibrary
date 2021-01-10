#pragma once
#include "vector_space.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include "complex.h"

namespace math_rz
{
	template<typename F,int n,int m>
	class matrix;
	template<typename F, int n>
	class finite_dimensional_vector_space :public vector_space<F>
	{
	public:

		finite_dimensional_vector_space(const std::vector<F>& a) :u(n)
		{
			if (n != a.size())
				throw std::domain_error("Dimensions are not compatible");
			std::copy(a.begin(), a.end(), u.begin());
		}
		finite_dimensional_vector_space(std::vector<F>&& a) :u(std::move(a))
		{
			if (n != u.size())
				throw std::domain_error("Dimensions are not compatible");
		}

		finite_dimensional_vector_space() :u(n) {}
		inline constexpr static int dimension = n;
		static finite_dimensional_vector_space _0()
		{
			return finite_dimensional_vector_space();
		}
		finite_dimensional_vector_space& operator+=(const finite_dimensional_vector_space& o)
		{
			for (int i = 0; i < n; i++)
				u.at(i) += o.u.at(i);
			return *this;
		}

		finite_dimensional_vector_space& operator-=(const finite_dimensional_vector_space& o)
		{
			for (int i = 0; i < n; i++)
				u.at(i) -= o.u.at(i);
			return *this;
		}
		finite_dimensional_vector_space& operator*=(const F& k)
		{
			for (int i = 0; i < n; i++)
				u.at(i) *= k;
			return *this;
		}

		finite_dimensional_vector_space& operator/=(const F& k)
		{
			for (int i = 0; i < n; i++)
				u.at(i) /= k;
			return *this;
		}

		finite_dimensional_vector_space& operator*=(int k)
		{
			for (int i = 0; i < n; i++)
				u.at(i) *= k;
			return *this;
		}
	/*	finite_dimensional_vector_space& operator/=(int k)
		{
			for (int i = 0; i < n; i++)
				u.at(i) /= k;
			return *this;
		}*/

		finite_dimensional_vector_space operator-() const
		{
			finite_dimensional_vector_space p;
			std::transform(u.begin(), u.end(), p.u.begin(), [](auto a) {return -a; });
			return p;
		}

		finite_dimensional_vector_space conj() const requires field_constraints::is_complex<F>
		{
			finite_dimensional_vector_space w = (*this);
			for (auto& s : w.u)
				s = s.conj();
			return w;
		}
		matrix<F, 1, n> transpose() const
		{
			return matrix<F, 1, n>({ this->u });
		}

		matrix<F, n, 1> as_matrix() const
		{
			return transpose().transpose();
		}
		template<int m>
		matrix<F, n, m> outer_product(const finite_dimensional_vector_space<F, m>&s) const
		{
			return as_matrix()*s.transpose();
		}

		template<int m>
		finite_dimensional_vector_space<F, n*m> kroenecker_product(const finite_dimensional_vector_space<F, m>& s) const
		{
			return (as_matrix() * s.transpose()).as_vector();
		}

		const F& operator[](int i) const
		{
			return u[i];
		}
		F& operator[](int i)
		{
			return u[i];
		}
		const F& at(int i) const
		{
			return u.at(i);
		}
		F& at(int i)
		{
			return u.at(i);
		}
		bool is_zero() const
		{
			return all_of(u.begin(), u.end(), [](const auto& x) {return x.is_zero(); });
		}

		std::vector<F>& get_vect()
		{
			return u;
		}

		const std::vector<F>& get_vect() const
		{
			return u;
		}
	protected:
		std::vector<F> u;
	};
	template <typename F, int n>
	using coordinate_space = finite_dimensional_vector_space<F, n>;

	template <typename F, int n>
	finite_dimensional_vector_space<F, n> operator+(
		const finite_dimensional_vector_space<F, n>& a, const finite_dimensional_vector_space<F, n>& b)
	{
		auto c(a);
		return c += b;
	}

	template <typename F, int n>
	finite_dimensional_vector_space<F, n> operator-(
		const finite_dimensional_vector_space<F, n>& a, const finite_dimensional_vector_space<F, n>& b)
	{
		auto c(a);
		return c -= b;
	}

	template <typename F, int n>
	finite_dimensional_vector_space<F, n> operator*(
		const F& k, const finite_dimensional_vector_space<F, n>& a)
	{
		auto c(a);
		return c *= k;
	}

	template <typename F, int n>
	finite_dimensional_vector_space<F, n> operator/(
		const F& k, const finite_dimensional_vector_space<F, n>& a)
	{
		auto c(a);
		return c /= k;
	}

	template <typename F, int n>
	finite_dimensional_vector_space<F, n> operator*(int k, const finite_dimensional_vector_space<F, n>& a)
	{
		auto c(a);
		return c *= k;
	}

	template <typename F, int n>
	finite_dimensional_vector_space<F, n> operator/(int k, const finite_dimensional_vector_space<F, n>& a)
	{
		auto c(a);
		return c /= k;
	}

	template <int n>
	finite_dimensional_vector_space<complex, n> operator*(
		const real_field& k, const finite_dimensional_vector_space<complex, n>& a)
	{
		auto c(a);
		return c *= {k, 0};
	}

	template <int n>
	finite_dimensional_vector_space<complex, n> operator/(
		const real_field& k, const finite_dimensional_vector_space<complex, n>& a)
	{
		auto c(a);
		return c /= {k, 0};
	}


	template <typename F, int n>
	std::ostream& operator<<(std::ostream& H, const finite_dimensional_vector_space<F, n>& p)
	{
		H << "( ";
		for (int i = 0; i < n; i++)
			if (i == n - 1) H << p.at(i) << " )";
			else H << p.at(i) << ", ";
		return H;
	}

	namespace vector_space_constraint
	{
		template<typename F, int n, typename M>
		concept is_vector = std::is_base_of_v<finite_dimensional_vector_space<F, n>, M>;
	}
}