#pragma once
#include "matrix.h"
#include "poly/polynomial.h"
#include "absalg/rational_extension.h"

namespace math_rz
{
		template<typename K, int n>
	class square_matrix :virtual public matrix<K, n, n>, virtual  public ring
	{
	public:
		square_matrix() {}
		square_matrix(const std::vector<std::vector<K>>& M) :
			matrix<K, n, n>(M) {}
		square_matrix(std::vector<std::vector<K>>&& M) :matrix<K, n, n>(std::move(M)) {}
		square_matrix(matrix<K, n, n>&& M) :matrix<K, n, n>(std::move(M))
		{	}
		square_matrix(const matrix<K, n, n>& M) :matrix<K, n, n>(M) {}
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

		square_matrix transpose() const
		{
			square_matrix T;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					T[j][i] = this->at(i).at(j);
			return T;
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

		square_matrix& operator*=(const K& k)
		{
			this->matrix<K, n, n>::operator*=(k);
			return *this;
		}

		square_matrix& operator/=(const K& k)
		{
			this->matrix<K, n, n>::operator/=(k);
			return *this;
		}

		static square_matrix _0()
		{
			return std::move(square_matrix());
		}
		static square_matrix _1()
		{
			return std::move(square_matrix(1));
		}
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
		K det() const
		{
			square_matrix M(std::move(this->row_echelon_form()));
			K d = 1;
			for (int i = 0; i < n; i++)
				if (M.at(i).at(i).is_zero())
					return 0;
				else d *= M.at(i).at(i);
			return d;
		}

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

		polynomial<K> caracteristic_polynomial() const
		{
			square_matrix<rational_extension<polynomial<K>>, n> J;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					if (i != j)
						J.at(i).at(j) = polynomial<K>(this->at(i).at(j));
					else J.at(i).at(j) = polynomial<K>({ this->at(i).at(j),-1 });
				}
			return (polynomial<K>)J.det();
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