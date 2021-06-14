#pragma once
#include "vector_space.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include "complex.h"
#include "structure/vector/inner_product.h"
#include <concepts>
#include <utility>

namespace math_rz::linalg
{
	namespace vector_space_constraint
	{
		template<typename K, int n, typename M>
		concept is_vector = std::is_base_of_v<finite_dimensional_vector_space<K, n>, M>;

		template<typename X>
		concept metric_space = requires(const X&u,const X&v)
		{
			u.distance(v);
		};
		
		template<typename M>
		concept vector_space = requires
		{
			typename M::base_field;
			M::dimension;
		};


		template<typename E>
		concept normed_vector_space = vector_space<E> && requires(const E & u)
		{
			u.norm();
		};

		template<typename H>
		concept inner_product_space = normed_vector_space<H> && requires(const H & u, const H & v)
		{
			u.inner_product(v);
		};


		template<typename E, typename F>
		concept vector_space_over_same_base_field =
			vector_space<E> && vector_space<F> &&
			std::is_same<typename E::base_field, typename F::base_field>::value;

		template<typename E, typename F>
		concept isomorpthic_vector_spaces = vector_space_over_same_base_field<E, F> && (E::dimension == F::dimension);

		template<vector_space E, vector_space F> requires vector_space_over_same_base_field<E, F>
			using product_space = finite_dimensional_vector_space<typename E::base_field, E::dimension + F::dimension>;
	}
	template<typename K,int n,int m>
	class matrix;
	template<typename K, int n>
	class square_matrix;

	template<typename K, int n>
	class finite_dimensional_vector_space :public vector_space<K>
	{
	protected:
		using structure_type = math_rz::linalg::structure::vector::metric_topology<K,n>;
	public:
		template<typename H>
		finite_dimensional_vector_space(const std::vector<H>& a) :u(n)
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

		finite_dimensional_vector_space(const finite_dimensional_vector_space& a) :u(a.u)
		{
		}

		finite_dimensional_vector_space() :u(n) {}
		inline constexpr static int dimension = n;
		using base_field = K;
		/*
		* Construct a vector from two vectors
		*/
		template<vector_space_constraint::vector_space E1, vector_space_constraint::vector_space E2>
		explicit finite_dimensional_vector_space(const E1& u1, const E2& u2) requires vector_space_constraint::vector_space_over_same_base_field<E1,E2>
			&& std::is_same_v<typename E1::base_field,K> && (E1::dimension+ E2::dimension ==n)
		{
			if constexpr (E1::dimension > 1) for (const auto& w : u1.get_vect())
				u.push_back(w);
			else u.push_back(u1);
			if constexpr (E2::dimension > 1) for (const auto& w : u2.get_vect())
				u.push_back(w);
			else u.push_back(u2);

		}

		/*
		* Construct a vector from three vectors
		*/
		template<vector_space_constraint::vector_space E1, vector_space_constraint::vector_space E2,
			vector_space_constraint::vector_space E3>
			explicit finite_dimensional_vector_space(const E1& u1, const E2& u2, const E3& u3)
			requires vector_space_constraint::vector_space_over_same_base_field<E1, E2>
			&& vector_space_constraint::vector_space_over_same_base_field<E2, E3>
			&& vector_space_constraint::vector_space_over_same_base_field<finite_dimensional_vector_space, E2> && (E1::dimension + E2::dimension + E3::dimension == n)
		{
			if constexpr (E1::dimension > 1) for (const auto& w : u1.get_vect())
				u.push_back(w);
			else u.push_back(u1);
			if constexpr (E2::dimension > 1) for (const auto& w : u2.get_vect())
				u.push_back(w);
			else u.push_back(u2);
			if constexpr (E3::dimension > 1) for (const auto& w : u3.get_vect())
				u.push_back(w);
			else u.push_back(u3);
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

		finite_dimensional_vector_space conj() const requires field_constraints::field_with_conj<K>
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
			return as_matrix().conj() *s.transpose();
		}

		template<int m>
		finite_dimensional_vector_space<K, n*m> kroenecker_product(const finite_dimensional_vector_space<K, m>& s) const
		{
			return (as_matrix() * s.transpose()).as_vector();
		}

		template<int p,int q> 
		std::conditional_t<p == q, square_matrix<K, p>, matrix<K, p, q>> reshape() const requires (p* q == n)
		{
			std::conditional_t<p == q, square_matrix<K, p>, matrix<K, p, q>>M;
			for (int i = 0; i < p; i++) for (int j = 0; j < q; j++)
				M[i][j] = u[i * q + j];
			return M;
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
		real_field metric(const finite_dimensional_vector_space& p) const
		{
			return structure_ptr->metric(*this, p);
		}

		real_field distance(const finite_dimensional_vector_space& p) const
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

		K dot_product(const finite_dimensional_vector_space& q) const
		{
			return dynamic_cast<math_rz::linalg::structure::vector::inner_product_topology<K, n>*>
				(structure_ptr.get())->dot_product(*this, q);
		}
	protected:
		std::vector<K> u;

		inline static std::unique_ptr<structure_type> structure_ptr =
			std::unique_ptr<structure_type>
			(new std::conditional_t<field_constraints::field_with_abs<K>,
				math_rz::linalg::structure::vector::L2_vect_inner_product<K, n>,
				math_rz::linalg::structure::vector::hamming_metric<K, n>>);
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
}