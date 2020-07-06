#pragma once
#include "matrix.h"
#include "polynomial.h"
#include "rational_extension.h"

template<typename F,int n>
class square_matrix :public matrix<F, n, n>, public ring
{
public:
	square_matrix(){}
	square_matrix(const std::vector<std::vector<F>>& M):
		matrix<F, n,n>(M){}
	square_matrix(std::vector<std::vector<F>>&& M):matrix<F, n, n>(std::move(M)) {}
	square_matrix(matrix<F, n, n>&& M) :matrix<F, n, n>(std::move(M)) 
	{	}
	square_matrix(const matrix<F, n, n>& M) :matrix<F, n, n>(M) {}
	square_matrix(int k):square_matrix()
	{
		for (int i = 0; i < n; i++)
			this->u[i][i] = k;

	}
	square_matrix& operator+=(const square_matrix& o)
	{
		return matrix<F, n, n>::operator+=(o);
	}
	square_matrix& operator-=(const square_matrix& o)
	{
		return matrix<F, n, n>::operator-=(o);
	}
	square_matrix operator-() const
	{
		return matrix<F, n, n>::operator-();
	}
	square_matrix& operator*=(const square_matrix& M)
	{
		square_matrix<F, n> P;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < n; k++)
					P.at(i).at(j) += this->at(i).at(k) * M.at(k).at(j);
		*this = std::move(P);
		return *this;
	}

	square_matrix& operator*=(const F& k)
	{
		return this->matrix<F,n,n>::operator*=(k);
	}

	square_matrix& operator/=(const F& k)
	{
		return this->matrix<F, n, n>::operator/=(k);
	}

	static square_matrix _0()
	{
		return std::move(square_matrix());
	}
	static square_matrix _1()
	{
		return std::move(square_matrix(1));
	}

	F det() const
	{
		square_matrix M (std::move(this->row_echelon_form()));
		F d = F::_1();
		for (int i = 0; i < n; i++)
			if (M.at(i).at(i) == F::_0())
				return F::_0();
			else d *= M.at(i).at(i);
		return d;
	}

	polynomial<F> caracteristic_polynomial() const
	{
		square_matrix<rational_extension<polynomial<F>>,n> J(square_matrix<rational_extension<polynomial<F>>,n>::_1());
		//for (int i = 0; i < n; i++)
			//for (int j = 0; j < n; j++)
				//J.at(i).at(j) -= polynomial<F>(this->at(i).at(j));
		return J.det().nominator();
	}
};

template <typename F, int n, int m>
square_matrix<F, n> operator+(
	const square_matrix<F, n>& a, const square_matrix<F, n>& b)
{
	auto c(a);
	return c += b;
}

template <typename F, int n, int m>
square_matrix<F, n> operator-(
	const square_matrix<F, n>& a, const square_matrix<F, n>& b)
{
	auto c(a);
	return c -= b;
}

template <typename F, int n, int m>
square_matrix<F, n> operator*(
	const F& k, const square_matrix<F,n>& a)
{
	auto c(a);
	return c *= k;
}

template <typename F, int n>
square_matrix <F, n> operator*(
	const square_matrix<F, n>& M, const square_matrix<F, n>& N)
{
	square_matrix<F, n> P;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				P.at(i).at(j) += M.at(i).at(k) * N.at(k).at(j);
	return P;
}