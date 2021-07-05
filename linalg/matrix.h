#pragma once
#include "finite_dimensional_vector_space.h"
#include "linalg/structure/matrix/inner_product.h"
#include <type_traits>
#include "absalg/ring.h"
#include "poly/polynomial.h"


namespace math_rz::linalg
{
	template<typename K, int n, int m = n>
	class matrix : virtual public std::conditional_t<n==m,ring,group>
	{
	protected:
		using structure_type = math_rz::linalg::structure::matrix::metric_topology<K, n, m>;
		inline static constexpr struct empty_matrix_t {} empty_matrix;
		matrix(empty_matrix_t) {};
	public:
		matrix() :u(n)
		{
			for (auto& v : u)
				v.resize(m);
		}
		matrix(const std::vector<std::vector<K>>& M) :
			u(M)
		{
			for (const auto& v : u)
				if (v.size() != m)
					throw std::domain_error("Dimensions are not compatible");
		}
		matrix(std::vector<std::vector<K>>&& M) :
			u(std::move(M))
		{
			for (const auto& v : u)
				if (v.size() != m)
					throw std::domain_error("Dimensions are not compatible");
		}

		template<typename H>
		matrix(const matrix<H, n, m>& M) :
			u(n, std::vector<K>(m))
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					u[i][j] = H(M[i][j]);
		}

		matrix(const K& k) requires (n == m) : matrix<K, n, m>()
		{
			for (int i = 0; i < n; i++)
				this->u[i][i] = k;
		}

		matrix(int k) requires(n == m) : matrix<K, n, m>()
		{
			for (int i = 0; i < n; i++)
				this->u[i][i] = k;

		}

		constexpr static int dimension = n * m;
		constexpr static int domain_dimension = m;
		constexpr static int codomain_dimension = n;
		using base_field = K;

		matrix<K, m, n> transpose() const
		{
			matrix<K, m, n> T;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					T[j][i] = this->at(i).at(j);
			return T;
		}
		matrix<K, m, n> T() const
		{
			return transpose();
		}

		matrix<K, m, n> conj_transpose() const
		{
			matrix<K, m, n> T;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					T[j][i] = this->at(i).at(j).conj();
			return T;
		}
		matrix<K, m, n> H() const
		{
			return conj_transpose();
		}

		matrix<K, m, n> hermitian_transpose() const
		{
			return conj_transpose();
		}

		matrix conj() const
		{
			matrix T;
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


		/*
		* matrix multiplication by a matrix
		*/
		matrix& operator*=(const matrix& M) requires(n==m)
		{
			*this = std::move((*this)*M);
			return *this;
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

		K trace() const
		{
			K s;
			for (int i = 0; i < std::min(n, m); i++)
				s += this->u[i][i];
			return s;
		}
		K tr() const
		{
			return trace();
		}

		bool is_one() const
		{
			if constexpr (n != m)
				return false;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					if (i != j)
					{
						if (!this->u[i][j].is_zero())
							return false;
					}
					else if (!this->u[i][i].is_one())
						return false;
			return true;
		}


		K det() const requires(n == m)
		{
			using matrix_type_extension = std::conditional_t < field_constraints::field<K>,
				matrix, matrix < rational_extension<K>, n, n >>;
			using ring_extension = std::conditional_t < field_constraints::field<K>,
				K, rational_extension<K>>;
			static_assert(field_constraints::field<rational_extension<K>>);
			matrix_type_extension S(*this);
			matrix_type_extension M(S.row_echelon_form());
			ring_extension d(1);
			for (int i = 0; i < n; i++)
				if (M.at(i).at(i).is_zero())
					return 0;
				else d *= M.at(i).at(i);
			return (K)d;
		}

		/*
		* Calculate the inverse of the given matrix
		* the complexity of this method is O(n^3)
		*/
		matrix inv() const requires(n == m)
		{
			auto Q = matrix<K, n, 2 * n>();
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					Q[i][j] = this->u[i][j];
			for (int i = 0; i < n; i++)
				for (int j = n; j < 2 * n; j++)
					Q[i][j] = i == (j - n);
			auto S = Q.row_echelon_form();
			matrix M1, M2;
			for (int i = 0; i < n; i++)
			{
				auto  r = S[i][i];
				if (r.is_zero())
					continue;
				for (int j = 0; j < n; j++)
				{
					M1[i][j] = S[i][j] / r;
					M2[i][j] = S[i][j + n] / r;
				}
			}
			for (int i = n - 1; i >= 0; i--)
			{
				for (int j = i - 1; j >= 0; j--)
				{
					auto r = M1[j][i];
					for (int k = 0; k < n; k++)
					{
						M2[j][k] -= r * M2[i][k];
						M1[j][k] -= r * M1[i][k];
					}
				}
			}

			return M2;
		}

		matrix  operator/=(const matrix& M) requires(n == m)
		{
			return *this *= M.inv();
		}

		/*
		* calculate the characteristic polynomial of this matrix
		* Note that this method is of exponential complexity, this is because of the exponential
		* blow up of degrees of the intermediate polynomials in the calculations
		* This effect is presented because of Guassian-Elimination
		* There exist a polynomial time algorithm for this problem
		*/
		poly::polynomial<K> caracteristic_polynomial() const requires(n == m)
		{
			matrix<rational_extension<poly::polynomial<K>>, n, n> J;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					if (i != j)
						J.at(i).at(j) = poly::polynomial<K>(this->at(i).at(j));
					else J.at(i).at(j) = poly::polynomial<K>({ this->at(i).at(j),-1 });
				}
			return (poly::polynomial<K>)J.det();
		}

		bool is_zero() const override
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					if (!u[i][j].is_zero())
						return false;
			return true;
		}


		matrix& operator*=(const K& k)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					u.at(i).at(j) *= k;
			return *this;
		}

		matrix& operator/=(const K& k)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					u.at(i).at(j) /= k;
			return *this;
		}

		matrix& operator/=(int k)
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
				for (; (i + off) < m && M[i][i + off].is_zero(); off++, r++);
			return r;
		}

		coordinate_space<K, n* m> as_vector() const
		{
			coordinate_space<K, n* m> p;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					p[i * m + j] = this->u[i][j];
			return p;
		}

		coordinate_space<K,n*m> flatten() const
		{
			return as_vector();
		}

		std::vector<std::vector<K>>& get_vect_vect()
		{
			return u;
		}

		const std::vector<std::vector<K>>& get_vect_vect() const
		{
			return u;
		}

		std::vector<K> get_column(int k) const
		{
			std::vector<K> S;
			for (int i = 0; i < n; i++)
				S.emplace_back(u[i][k]);
			return S;
		}


		std::vector<std::reference_wrapper<K>> get_column_ref(int k)
		{
			std::vector<std::reference_wrapper<K>> S;
			for (int i = 0; i < n; i++)
				S.push_back(u[i][k]);
			return S;
		}

		static void set_structure(structure_type* S)
		{
			structure_ptr.reset(S);
		}
		static const structure_type& get_structure()
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
			return dynamic_cast<math_rz::linalg::structure::matrix::norm_topology<K, n, m>*>
				(structure_ptr.get())->norm(*this);
		}

		K inner_product(const matrix& q) const
		{
			return dynamic_cast<math_rz::linalg::structure::matrix::inner_product_topology<K, n, m>*>
				(structure_ptr.get())->inner_product(*this, q);
		}

		K dot_product(const matrix& q) const
		{
			return dynamic_cast<math_rz::linalg::structure::matrix::inner_product_topology<K, n, m>*>
				(structure_ptr.get())->dot_product(*this, q);
		}

		void foreach(const std::function<void(K&)>& f)
		{
			for (auto& R : u) for (auto& a : R)
				f(a);
		}

		finite_dimensional_vector_space<K, m> solve(finite_dimensional_vector_space<K, m> b) const
		{
			matrix<K, n, m> C;
			real_field sigma = 1;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					C[i][j] = u[i][j];
			for (int i = 0; i < n; i++)
			{
				int p;
				for (p = i; p < n && C[i][p] == 0; p++);
				if (p != i && p < n)
				{
					std::iter_swap(C.begin() + i, C.begin() + p);
					std::iter_swap(b.get_vect().begin() + i, b.get_vect().begin() + p);
				}
				else if (p == n)
					break;

				for (int j = i + 1; j < n; j++)
				{
					b[j] -= (C[j][i] / C[i][i]) * b[i];
					C[j] -= (C[j][i] / C[i][i]) * C[i];
				}
			}
			for (int i = n - 1; i >= 0; i--) if (C[i][i] != 0)
			{
				b[i] /= C[i][i];
				C[i] /= C[i][i];
				for (int j = 0; j < i; j++)
				{
					b[j] -= C[j][i] * b[i];
					C[j] -= C[j][i] * C[i];
				}
			}
			return b;
		}
		matrix<K, n, m> row_echelon_form(bool down = false) const
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
					K r = P[j][p] / P[i][p];
					for (int k = p; k < m; k++)
						P[j][k] = P[j][k] - r * P[i][k];
				}
			}
			return P;
		}
		operator K() const requires(n == m && n == 1)
		{
			return u[0][0];
		}

		operator finite_dimensional_vector_space<K, n>() const requires(n>1 && m == 1)
		{
			finite_dimensional_vector_space<K, n> S;
			for (int i = 0; i < n; i++)
				S[i] = u[i][0];
			return S;
		}
	protected:
		//boost::multi_array<K,2> u;
		//using structure_type = math_rz::linalg::structure::matrix::L22_operator_norm;
		std::vector<std::vector<K>> u;

		/*
		* This variable gives the analytical structure of the matrix
		* If the underlying field has a norm, use the L22 operator norm
		* else use the hamming metric
		*/
		inline static std::unique_ptr<structure_type> structure_ptr =
			std::unique_ptr<structure_type>
			(new std::conditional_t< field_constraints::field_with_abs<K>,
				math_rz::linalg::structure::matrix::L22_operator_norm<K, n, m>,
				math_rz::linalg::structure::matrix::hamming_metric<K, n, m>>);
	};

	namespace matrix_constraint {
		template<typename K, int n, int m, typename M>
		concept is_matrix = std::is_base_of_v<matrix<K, n, m>, M>;
	}

	template <typename K, int n, int m>
	matrix<K, n, m> operator+(
		const matrix<K, n, m>& a, const matrix<K, n, m>& b)
	{
		auto c(a);
		return c += b;
	}

	template <typename K, int n, int m>
	matrix<K, n, m> operator-(
		const matrix<K, n, m>& a, const matrix<K, n, m>& b)
	{
		auto c(a);
		return c -= b;
	}

	template <typename K, int n, int m>
	matrix<K, n, m> operator*(
		const K& k, const matrix<K, n, m>& a)
	{
		auto c(a);
		return c *= k;
	}

	template <typename K, int n, int p, int m>
	matrix<K, n, m> operator*(
		const matrix<K, n, p>& M, const matrix<K, p, m>& N)
	{
		/*
		* This nested for loop will calculate the matrix product
		* The order of the last two for is intentionally inverted to reduce cache misses
		*/
		matrix<K, n, m> P;
		for (int i = 0; i < n; i++)
			for (int k = 0; k < p; k++)
				for (int j = 0; j < m; j++)
					P.at(i).at(j) += M.at(i).at(k) * N.at(k).at(j);
		return P;
	}

	template <typename K, int n>
	K operator*(
		const matrix<K, 1, n>& M, const matrix<K, n, 1>& N)
	{
		K P;
		for (int i = 0; i < n; i++)
			P += M.at(0).at(i) * N.at(i).at(0);
		return P;
	}

	template <typename K, int n>
	K operator*(
		const matrix<K, 1, n>& M, const coordinate_space<K, n>& N)
	{
		K P;
		for (int i = 0; i < n; i++)
			P += M.at(0).at(i) * N.at(i);
		return P;
	}

	template <typename K, int n, int m>
	matrix<K, n, m> operator*(const matrix<K, n, m>& M, const K& k)
	{
		auto c(M);
		return c *= k;
	}

	template <typename K, int n, int m>
	coordinate_space<K, n> operator*(const matrix<K, n, m>& M, const coordinate_space<K, m>& u)
	{
		coordinate_space<K, n> v;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				v.at(i) += M.at(i).at(j) * u.at(j);
		return v;
	}

	template <typename K, int n, int m>
	matrix<K, n, m> operator/(const matrix<K, n, m>& M, const K& k)
	{
		auto c(M);
		return c /= k;
	}

	template <typename K, int n, int m>
	matrix<K, n, m> operator/(const matrix<K, n, m>& M,int k)
	{
		auto c(M);
		return c /= k;
	}

	template <typename K, int n, int m>
	std::ostream& operator<<(std::ostream& H, const matrix<K, n, m>& p)
	{

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
				H << p.at(i).at(j) << '\t';
			H << '\n';
		}
		return H;
	}

	template <typename K, int n, int m>
	std::istream& operator>>(std::istream& H, matrix<K, n, m>& p)
	{

		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				H >> p.at(i).at(j);
		return H;
	}
}