#pragma once
#include "finite_dimensional_vector_space.h"

template<typename F,int n,int m>
class matrix :public finite_dimensional_vector_space<std::vector<F> , n >
{
public:
	matrix()
	{
		for (auto& v : u)
			v.resize(m);
	}
	matrix(const std::vector<std::vector<F>> &M):
		finite_dimensional_vector_space<F,n>(M)
	{
		for (const auto& v : u)
			if (v.size() != m)
				throw std::domain_error("Dimensions are not compatible");
	}
	matrix(std::vector<std::vector<F>>&& M) :
		finite_dimensional_vector_space<std::vector<F>, n>(std::move(M))
	{
		for (const auto& v : u)
			if (v.size() != m)
				throw std::domain_error("Dimensions are not compatible");
	}

	matrix<F, m, n> transpose() const
	{
		matrix<F, m, n> T;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				T[j][i] = this->at(i).at(j);
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
			for(int j=0;j<m;j++)
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
	matrix row_echelon_form(bool down=false) const
	{
		matrix P = (*this);
		for (int i = 0, p = 0; (i < n)&&(p<m); i++,p++)
		{
			if (P.at(i).at(p).is_zero())
			{
				int j;
				for (j = i + 1; j < n; j++)
					if (P.at(j).at(p) != F::_0())
					{
						std::swap(P.at(i), P.at(j));
						break;
					}
				if (j == n)
					continue;
				
			}
			for (int j = i + 1; j < n; j++)
			{
					F&& r = P.at(j).at(p) / P.at(i).at(p);
				for (int k = p; k < m; k++)
					P.at(j).at(k) = P.at(j).at(k) - r * P.at(i).at(k);
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
		int r=0;
		for (int i = n - 1; i >= 0; i--)
			if (any_of(u.at(i).rbegin(), u.at(i).rend(), [](const auto& c)
				{
					return c != 0;
				}))
				return r;
			else r++;
		return r;
	}
	bool is_zero() const
	{
		for(const auto& v:u)
			if (any_of(v.begin(), v.end(), [](auto& x) 
				{
					return !x.is_zero();
				}))
				return false;
		return true;
	}
protected:
	using finite_dimensional_vector_space<std::vector<F>, n > ::u;
};

template <typename F, int n,int m>
matrix<F, n,m> operator+(
	const matrix<F, n,m>& a, const matrix<F, n,m>& b)
{
	auto c(a);
	return c += b;
}

template <typename F, int n, int m>
matrix<F, n,m> operator-(
	const matrix<F, n,m>& a, const matrix<F, n,m>& b)
{
	auto c(a);
	return c -= b;
}

template <typename F, int n, int m>
matrix<F, n,m> operator*(
	const F& k, const matrix<F, n,m>& a)
{
	auto c(a);
	return c *= k;
}

template <typename F, int n,int p, int m>
matrix<F, n, m> operator*(
	const matrix<F, n, p>& M, const matrix<F, p, m>& N)
{
	matrix<F,n,m> P;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			for (int k = 0; k < p; k++)
				P.at(i).at(j) += M.at(i).at(k) * N.at(k).at(j);
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
coordinate_space<F, n> operator*(const matrix<F, n, m>& M,const coordinate_space<F,m> &u)
{
	coordinate_space<F,n> v;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			v.at(i) += M.at(i).at(j) * u.at(j);
	return v;
}

template <typename F, int n, int m>
matrix<F, n, m> operator/(const matrix<F, n, m>& M,const F& k)
{
	auto c(M);
	return c /= k;
}

template <typename F, int n,int m>
std::ostream& operator<<(std::ostream& H, const matrix<F, n,m>& p)
{
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			H << p.at(i).at(j) << '\t';
		H << '\n';
	}
	return H;
}