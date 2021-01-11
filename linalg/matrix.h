#pragma once
#include "finite_dimensional_vector_space.h"
#include "boost/multi_array.hpp"
#include "linalg/structure/matrix/inner_product.h"


namespace math_rz {
	template<typename F, int n, int m>
	class matrix: virtual public group
	{
	protected:
		using structure_type = math_rz::linalg::structure::matrix::metric_topology<F, n, m>;
	public:
		matrix():u(n)
		{
			for (auto& v : u)
				v.resize(m);
		}
		matrix(const std::vector<std::vector<F>>& M) :
			u(M)
		{
			for (const auto& v : u)
				if (v.size() != m)
					throw std::domain_error("Dimensions are not compatible");
		}
		matrix(std::vector<std::vector<F>>&& M) :
			u(std::move(M))
		{
			for (const auto& v : u)
				if (v.size() != m)
					throw std::domain_error("Dimensions are not compatible");
		}
		constexpr static int dimension = n * m;
		constexpr static int domain_dimension = m;
		constexpr static int codomain_dimension = n;
		using base_field = F;

		matrix<F, m, n> transpose() const
		{
			matrix<F, m, n> T;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					T[j][i] = this->at(i).at(j);
			return T;
		}

		matrix<F, m, n> conj_transpose() const
		{
			matrix<F, m, n> T;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					T[j][i] = this->at(i).at(j).conj();
			return T;
		}

		matrix<F, m, n> conj() const
		{
			matrix<F, m, n> T;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					T[i][j] = this->at(i).at(j).conj();
			return T;
		}

		bool operator!=(const matrix& M) const
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					if (this->at(i).at(j) != M.at(i).at(j))
						return true;
			return false;
		}

		bool operator==(const matrix& M) const
		{
			return !(*this != M);
		}
		matrix& operator+=(const matrix& o)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					u.at(i).at(j) += o.u.at(i).at(j);
			return *this;
		}
		matrix& operator-=(const matrix& o)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					u.at(i).at(j) -= o.u.at(i).at(j);
			return *this;
		}
		matrix operator-() const
		{
			matrix p;
			std::transform(u.begin(), u.end(), p.u.begin(), [](auto a)
				{
					std::transform(a.begin(), a.end(), a.begin(), [](auto b)
						{
							return -b;
						});
					return a;
				});
			return p;
		}

		F trace() const
		{
			F tr;
			for (int i = 0; i < std::min(n, m); i++)
				tr += this->u[i][i];
			return tr;
		}
		F tr() const
		{
			return trace();
		}

		matrix& operator*=(const F& k)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					u.at(i).at(j) *= k;
			return *this;
		}

		matrix& operator/=(const F& k)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					u.at(i).at(j) /= k;
			return *this;
		}

		auto& operator[](int i)
		{
			return u[i];
		}
		const auto& operator[](int i) const
		{
			return u[i];
		}

		auto& at(int i)
		{
			return u.at(i);
		}
		const auto& at(int i) const
		{
			return u[i];
		}


		matrix row_echelon_form(bool down = false) const
		{
			matrix P = (*this);
			for (int i = 0, p = 0; (i < n) && (p < m); i++, p++)
			{
				if (P[i][p].is_zero())
				{
					int j;
					for (j = i + 1; j < n; j++)
						if (!P[j][p].is_zero())
						{
							std::swap(P[i], P[j]);
							break;
						}
					if (j == n)
						continue;

				}
				for (int j = i + 1; j < n; j++)
				{
					F&& r = P[j][p] / P[i][p];
					for (int k = p; k < m; k++)
						P[j][k] = P[j][k] - r * P[i][k];
				}
			}
			return P;
		}
		int rank() const
		{
			return m - nullity();
		}
		int nullity() const
		{
			int r = 0;
			auto M = row_echelon_form();
			int off = 0;
			for (int i = 0; i < std::min(n, m); i++)
				for (;(i+off)<m &&M[i][i + off].is_zero(); off++, r++);
			return r;
		}
		bool is_zero() const override
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					if (!u[i][j].is_zero())
						return false;
			return true;
		}

		coordinate_space<F, n* m> as_vector() const
		{
			coordinate_space<F, n* m> p;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					p[i * m + j] = this->u[i][j];
			return p;
		}

		std::vector<std::vector<F>>& get_vect_vect()
		{
			return u;
		}

		const std::vector<std::vector<F>>& get_vect_vect() const
		{
			return u;
		}
		static void set_structure(structure_type* S)
		{
			structure_ptr.reset(S);
		}
		static const structure_type & get_structure()
		{
			return (*structure_ptr);
		}

		real_field metric(const matrix& p)
		{
			return structure_ptr->metric(*this, p);
		}

		real_field distance(const matrix& p)
		{
			return metric(p);
		}

		real_field norm() const
		{
			return dynamic_cast<math_rz::linalg::structure::matrix::norm_topology<F,n,m>*>
				(structure_ptr.get())->norm(*this);
		}

		F inner_product(const matrix& q) const
		{
			return dynamic_cast<math_rz::linalg::structure::matrix::inner_product_topology<F,n,m>*>
				(structure_ptr.get())->inner_product(*this, q);
		}
	protected:
		//boost::multi_array<F,2> u;
		//using structure_type = math_rz::linalg::structure::matrix::L22_operator_norm;
		std::vector<std::vector<F>> u;
		inline static std::unique_ptr<structure_type> structure_ptr =
			std::unique_ptr<structure_type>
			(new math_rz::linalg::structure::matrix::L22_operator_norm<F,n,m>);
	};

	namespace matrix_constraint {
		template<typename F, int n, int m, typename M>
		concept is_matrix = std::is_base_of_v<matrix<F, n, m>, M>;
	}

	template <typename F, int n, int m>
	matrix<F, n, m> operator+(
		const matrix<F, n, m>& a, const matrix<F, n, m>& b)
	{
		auto c(a);
		return c += b;
	}

	template <typename F, int n, int m>
	matrix<F, n, m> operator-(
		const matrix<F, n, m>& a, const matrix<F, n, m>& b)
	{
		auto c(a);
		return c -= b;
	}

	template <typename F, int n, int m>
	matrix<F, n, m> operator*(
		const F& k, const matrix<F, n, m>& a)
	{
		auto c(a);
		return c *= k;
	}

	template <typename F, int n, int p, int m>
	matrix<F, n, m> operator*(
		const matrix<F, n, p>& M, const matrix<F, p, m>& N)
	{
		/*
		* This nested for loop will calculate the matrix product
		* The order of the last two for is intentionally inverted to reduce cache misses
		*/
		matrix<F, n, m> P;
		for (int i = 0; i < n; i++)
			for (int k = 0; k < p; k++)
				for (int j = 0; j < m; j++)
					P.at(i).at(j) += M.at(i).at(k) * N.at(k).at(j);
		return P;
	}

	template <typename F, int n>
	F operator*(
		const matrix<F, 1,n>& M, const matrix<F, n, 1>& N)
	{
		F P;
		for (int i = 0; i < n; i++)
			P += M.at(0).at(i) * N.at(i).at(0);
		return P;
	}

	template <typename F, int n>
	F operator*(
		const matrix<F, 1, n>& M, const coordinate_space<F,n>& N)
	{
		F P;
		for (int i = 0; i < n; i++)
			P += M.at(0).at(i) * N.at(i);
		return P;
	}

	template <typename F, int n, int m>
	matrix<F, n, m> operator/(const F& k, const matrix<F, n, m>& M)
	{
		auto c(M);
		return c /= k;
	}

	template <typename F, int n, int m>
	matrix<F, n, m> operator*(const matrix<F, n, m>& M, const F& k)
	{
		auto c(M);
		return c *= k;
	}

	template <typename F, int n, int m>
	coordinate_space<F, n> operator*(const matrix<F, n, m>& M, const coordinate_space<F, m>& u)
	{
		coordinate_space<F, n> v;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				v.at(i) += M.at(i).at(j) * u.at(j);
		return v;
	}

	template <typename F, int n, int m>
	matrix<F, n, m> operator/(const matrix<F, n, m>& M, const F& k)
	{
		auto c(M);
		return c /= k;
	}

	template <typename F, int n, int m>
	std::ostream& operator<<(std::ostream& H, const matrix<F, n, m>& p)
	{

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				H << p.at(i).at(j) << '\t';
			H << '\n';
		}
		return H;
	}

	template <typename F, int n, int m>
	std::istream& operator>>(std::istream& H, matrix<F, n, m>& p)
	{

		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				H >> p.at(i).at(j);
		return H;
	}
}