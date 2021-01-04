#pragma once
#include "derivator.h"

namespace math_rz
{

	/*
	* General Case:
	*/
	template<typename F, int m,int n,
		typename E1,
		typename E2>
		class default_derivator : public derivator<F, m, n, E1, E2>
	{
	public:
		default_derivator(E1 p0, F _eps) :derivator<F, m, n, E1, E2>(p0), eps(_eps) {}

		matrix<F, n, m> jacobian(const function<E1, E2>& f) const override
		{
			E1 s = this->x0;
			matrix<F, n, m> M;
			for (int i = 0; i < m; i++)
			{
				s[i] += eps;
				F k = F(1) / eps;
				E2 h = k * (f(s) - f(this->x0));
				s[i] -= eps;
				for (int j = 0; j < n; j++)
					M[j][i] = h[j];


			}
			return M;
		}

	private:
		F eps;
	};

	/*
	* When the two spaces are isomorphic (have the same dimension)
	*/

	template<typename F, int n,
		typename E1,
		typename E2>
		class default_derivator<F, n, n, E1, E2> : public derivator<F, n, n, E1, E2>
	{
	public:
		default_derivator(E1 p0, F _eps) :derivator<F, n, n, E1, E2>(p0), eps(_eps) {}

		square_matrix<F, n> jacobian(const function<E1, E2>& f) const override
		{
			E1 s = this->x0;
			square_matrix<F, n> M;
			for (int i = 0; i < n; i++)
			{
				s[i] += eps;
				F k = F(1) / eps;
				E2 h = k * (f(s) - f(this->x0));
				s[i] -= eps;
				for (int j = 0; j < n; j++)
					M[j][i] = h[j];


			}
			return M;
		}

	private:
		F eps;
	};

	/*
	* When the domain is the base field
	*/
	template<typename F, int n,
		typename E2>
		class default_derivator<F, 1, n, F, E2> : public derivator<F, 1, n, F, E2>
	{
	public:
		default_derivator(F p0, F _eps) :derivator<F, 1,n, F, E2>(p0), eps(_eps) {}

		matrix<F, n, 1> jacobian(const function<F, E2>& f) const override
		{
			F s = this->x0;
			matrix < F, n,1 > M;
			s += eps;
			F k = F(1) / eps;
			E2 h = k * (f(s) - f(this->x0));
			for (int j = 0; j < n; j++)
				M[j][0] = h[j];
			return M;
		}

	private:
		F eps;
	};


	/*
	* When the co-domain is the field F
	*/
	template<typename F, int m,
		typename E1>
		class default_derivator < F, m, 1, E1,F> : public derivator<F, m, 1, E1, F>
	{
	public:
		default_derivator(E1 p0, F _eps) :derivator<F, m, 1, E1, F>(p0), eps(_eps) {}

		matrix<F, 1, m> jacobian(const function<E1, F>& f) const override
		{
			E1 s = this->x0;
			matrix<F, 1, m> M;
			for (int i = 0; i < m; i++)
			{
				s[i] += eps;
				F k = F(1) / eps;
				F h = k * (f(s) - f(this->x0));
				s[i] -= eps;
				M[0][i] = h;


			}
			return M;
		}

	private:
		F eps;
	};

	/*
	* When both the domain and the co-domain are one dimensional
	*/

	template<typename F>
		class default_derivator<F, 1, 1, F,F> : public derivator<F, 1, 1, F, F>
	{
	public:
		default_derivator(F p0, F _eps) :derivator<F, 1, 1, F, F>(p0), eps(_eps) {}

		square_matrix<F, 1> jacobian(const function<F, F>& f) const override
		{
			F s = this->x0;
			matrix < F, 1, 1 > M;
			s += eps;
			F k = F(1) / eps;
			F h = k * (f(s) - f(this->x0));
			M[0][0] = h;
			return M;
		}

	private:
		F eps;
	};

}