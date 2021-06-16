#pragma once
#include "matrix.h"
#include "poly/polynomial.h"
#include "absalg/rational_extension.h"
#include <type_traits>

namespace math_rz::linalg
{
		template<typename K, int n>
	class square_matrix :virtual public matrix<K, n, n>, virtual  public ring
	{
	protected:
		square_matrix(matrix<K, n, n>::empty_matrix_t A):matrix<K,n,n>(A) {}
	public:
		square_matrix(const std::vector<K>& D)
		{
			if (D.size() != n)
				throw std::exception("Size Mismatch");
			for (int i = 0; i < n; i++)
				this->u[i][i] = D[i];
		}

		square_matrix(const coordinate_space<K,n>& D)
		{
			for (int i = 0; i < n; i++)
				this->u[i][i] = D[i];
		}
		square_matrix() {}
		square_matrix(const std::vector<std::vector<K>>& M) :
			matrix<K, n, n>(M) {}
		square_matrix(std::vector<std::vector<K>>&& M) :matrix<K, n, n>(std::move(M)) {}
		square_matrix(matrix<K, n, n>&& M) :matrix<K, n, n>(std::move(M))
		{	}
		square_matrix(const matrix<K, n, n>& M) :matrix<K, n, n>(M) {}
		template<typename H>
		square_matrix(const matrix<H,n,n>&M):matrix<K,n,n>(M){}
		square_matrix(const K &k) :square_matrix()
		{
			for (int i = 0; i < n; i++)
				this->u[i][i] = k;
		}

		square_matrix(int k) :square_matrix()
		{
			for (int i = 0; i < n; i++)
				this->u[i][i] = k;

		}

		template<typename H>
		square_matrix(const square_matrix<H, n>& M)
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					this->u[i][j] = M[i][j];
		}

		/*
		* returns the transpose of the matrix
		*/
		square_matrix transpose() const
		{
			square_matrix T;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					T[j][i] = this->at(i).at(j);
			return T;
		}

		/*
		* returns the Transpose of a matrix
		*/
		square_matrix T() const
		{
			return transpose();
		}

		/*
		* returns the Hermitian of a matrix
		*/
		square_matrix conj_transpose() const
		{
			square_matrix M(std::move(matrix<K, n, n>::conj_transpose()));
			return M;
		}

		/*
		* returns the Hermitian
		*/
		square_matrix H() const
		{
			return conj_transpose();
		}

		/*
		* returns the conjugate
		*/
		square_matrix conj() const
		{
			square_matrix M(std::move(matrix<K,n,n>::conj()));
			return M;
		}

		square_matrix& operator+=(const square_matrix& o)
		{
			matrix<K, n, n>::operator+=(o);
			return *this;
		}

		square_matrix& operator-=(const square_matrix& o)
		{
			matrix<K, n, n>::operator-=(o);
			return *this;
		}

		square_matrix operator-() const
		{
			matrix<K, n, n>::operator-();
			return *this;
		}

		/*
		* matrix multiplication by a matrix
		*/
		square_matrix& operator*=(const square_matrix& M)
		{
			square_matrix<K, n> P;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					for (int k = 0; k < n; k++)
						P.at(i).at(j) += this->at(i).at(k) * M.at(k).at(j);
			*this = std::move(P);
			return *this;
		}

		/*
		* Matrix multiplication by a scalar
		*/
		square_matrix& operator*=(const K& k)
		{
			this->matrix<K, n, n>::operator*=(k);
			return *this;
		}

		/*
		* matrix division by a scalar
		*/
		square_matrix& operator/=(const K& k)
		{
			this->matrix<K, n, n>::operator/=(k);
			return *this;
		}

		/*
		* is equal to the identity matrix
		*/
		bool is_one() const
		{
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


		/*
		* Calculates the determinant of the given matrix
		* The complexity of this method is O(n^3)
		*/
		K det() const
		{
			using matrix_type_extension = std::conditional_t<field_constraints::field<K>,
				square_matrix, square_matrix<rational_extension<K>,n>>;
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
		square_matrix inv() const
		{
			auto Q = matrix<K, n, 2 * n>();
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					Q[i][j] = this->u[i][j];
			for (int i = 0; i < n; i++)
				for (int j = n; j < 2 * n; j++)
					Q[i][j] = i == (j - n);
			auto S = Q.row_echelon_form();
			square_matrix M1, M2;
			for (int i = 0; i < n; i++)
			{
				auto  r = S[i][i];
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

		square_matrix operator/=(const square_matrix& M)
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
		poly::polynomial<K> caracteristic_polynomial() const
		{
			square_matrix<rational_extension<poly::polynomial<K>>, n> J;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					if (i != j)
						J.at(i).at(j) = poly::polynomial<K>(this->at(i).at(j));
					else J.at(i).at(j) = poly::polynomial<K>({ this->at(i).at(j),-1 });
				}
			return (poly::polynomial<K>)J.det();
		}
	};

	template <typename K, int n>
	square_matrix<K, n> operator+(
		const square_matrix<K, n>& a, const square_matrix<K, n>& b)
	{
		auto c(a);
		return c += b;
	}

	template <typename K, int n>
	square_matrix<K, n> operator-(
		const square_matrix<K, n>& a, const square_matrix<K, n>& b)
	{
		auto c(a);
		return c -= b;
	}

	template <typename K, int n>
	square_matrix <K, n> operator*(
		const square_matrix<K, n>& M, const square_matrix<K, n>& N)
	{
		square_matrix<K, n> P;
		for (int i = 0; i < n; i++)
			for (int k = 0; k < n; k++)
				for (int j = 0; j < n; j++)
					P.at(i).at(j) += M.at(i).at(k) * N.at(k).at(j);
		return P;
	}

	template <typename K, int n>
	square_matrix <K, n> operator/(
		const square_matrix<K, n>& M, const square_matrix<K, n>& N)
	{
		return M * N.inv();
	}


	template <typename K, int n>
	square_matrix <K, n> operator*(
		const K& k, const square_matrix<K, n>& N)
	{
		square_matrix<K, n> P;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				P.at(i).at(j) = k* N.at(i).at(j);
		return P;
	}
}