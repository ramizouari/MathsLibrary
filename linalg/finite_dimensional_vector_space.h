#pragma once
#include "vector_space.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include "complex.h"
#include "structure/vector/inner_product.h"

namespace math_rz
{
	template<typename K,int n,int m>
	class matrix;
	template<typename K, int n>
	class finite_dimensional_vector_space :public vector_space<K>
	{
	protected:
		using structure_type = math_rz::linalg::structure::vector::metric_topology<K,n>;
	public:

		finite_dimensional_vector_space(const std::vector<K>& a) :u(n)
		{
			if (n != a.size())
				throw std::domain_error("Dimensions are not compatible");
			std::copy(a.begin(), a.end(), u.begin());
		}
		finite_dimensional_vector_space(std::vector<K>&& a) :u(std::move(a))
		{
			if (n != u.size())
				throw std::domain_error("Dimensions are not compatible");
		}

		finite_dimensional_vector_space() :u(n) {}
		inline constexpr static int dimension = n;
		using base_field = K;
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
		finite_dimensional_vector_space& operator*=(const K& k)
		{
			for (int i = 0; i < n; i++)
				u.at(i) *= k;
			return *this;
		}

		finite_dimensional_vector_space& operator/=(const K& k)
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

		finite_dimensional_vector_space conj() const requires field_constraints::is_complex<K>
		{
			finite_dimensional_vector_space w = (*this);
			for (auto& s : w.u)
				s = s.conj();
			return w;
		}
		matrix<K, 1, n> transpose() const
		{
			return matrix<K, 1, n>({ this->u });
		}

		matrix<K, n, 1> as_matrix() const
		{
			return transpose().transpose();
		}
		template<int m>
		matrix<K, n, m> outer_product(const finite_dimensional_vector_space<K, m>&s) const
		{
			return as_matrix()*s.transpose();
		}

		template<int m>
		finite_dimensional_vector_space<K, n*m> kroenecker_product(const finite_dimensional_vector_space<K, m>& s) const
		{
			return (as_matrix() * s.transpose()).as_vector();
		}

		const K& operator[](int i) const
		{
			return u[i];
		}
		K& operator[](int i)
		{
			return u[i];
		}
		const K& at(int i) const
		{
			return u.at(i);
		}
		K& at(int i)
		{
			return u.at(i);
		}
		bool is_zero() const
		{
			return all_of(u.begin(), u.end(), [](const auto& x) {return x.is_zero(); });
		}

		std::vector<K>& get_vect()
		{
			return u;
		}

		const std::vector<K>& get_vect() const
		{
			return u;
		}

		static void set_structure(structure_type* S)
		{
			structure_ptr.reset(S);
		}
		static const structure_type& get_structure()
		{
			return (*structure_ptr);
		}
		real_field metric(const finite_dimensional_vector_space& p)
		{
			return structure_ptr->metric(*this, p);
		}

		real_field distance(const finite_dimensional_vector_space& p)
		{
			return metric(p);
		}

		real_field norm() const
		{
			return dynamic_cast<math_rz::linalg::structure::vector::norm_topology<K, n>*>
				(structure_ptr.get())->norm(*this);
		}

		K inner_product(const finite_dimensional_vector_space& q) const
		{
			return dynamic_cast<math_rz::linalg::structure::vector::inner_product_topology<K, n>*>
				(structure_ptr.get())->inner_product(*this, q);
		}
	protected:
		std::vector<K> u;

		inline static std::unique_ptr<structure_type> structure_ptr =
			std::unique_ptr<structure_type>
			(new math_rz::linalg::structure::vector::L2_vect_inner_product<K, n>);
	};
	template <typename K, int n>
	using coordinate_space = finite_dimensional_vector_space<K, n>;

	template <typename K, int n>
	finite_dimensional_vector_space<K, n> operator+(
		const finite_dimensional_vector_space<K, n>& a, const finite_dimensional_vector_space<K, n>& b)
	{
		auto c(a);
		return c += b;
	}

	template <typename K, int n>
	finite_dimensional_vector_space<K, n> operator-(
		const finite_dimensional_vector_space<K, n>& a, const finite_dimensional_vector_space<K, n>& b)
	{
		auto c(a);
		return c -= b;
	}

	template <typename K, int n>
	finite_dimensional_vector_space<K, n> operator*(
		const K& k, const finite_dimensional_vector_space<K, n>& a)
	{
		auto c(a);
		return c *= k;
	}

	template <typename K, int n>
	finite_dimensional_vector_space<K, n> operator/(
		const K& k, const finite_dimensional_vector_space<K, n>& a)
	{
		auto c(a);
		return c /= k;
	}

	template <typename K, int n>
	finite_dimensional_vector_space<K, n> operator*(int k, const finite_dimensional_vector_space<K, n>& a)
	{
		auto c(a);
		return c *= k;
	}

	template <typename K, int n>
	finite_dimensional_vector_space<K, n> operator/(int k, const finite_dimensional_vector_space<K, n>& a)
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


	template <typename K, int n>
	std::ostream& operator<<(std::ostream& H, const finite_dimensional_vector_space<K, n>& p)
	{
		H << "( ";
		for (int i = 0; i < n; i++)
			if (i == n - 1) H << p.at(i) << " )";
			else H << p.at(i) << ", ";
		return H;
	}

	namespace vector_space_constraint
	{
		template<typename K, int n, typename M>
		concept is_vector = std::is_base_of_v<finite_dimensional_vector_space<K, n>, M>;

		template<typename M>
		concept vector_space = requires 
		{
			typename M::base_field;
			M::dimension;
		};

		template<typename E, typename K>
		concept vector_space_over_same_base_field = 
			vector_space<E> && vector_space<K>&&
			std::is_same<typename E::base_field,typename K::base_field>::value;

		template<typename E,typename K>
		concept isomorpthic_vector_spaces = vector_space_over_same_base_field<E, K> && (E::dimension == K::dimension);

	}
}