#pragma once
#include "linalg/matrix.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "analysis/function.h"

/*
* - Let F be a field with an absolute value
* - Let E1 be an m dimensional vector space over F
* - Let E2 be an n dimensional vector space over F
* - We will consider the differentiable functions from E1 to E2 
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* This file contains the derivator template class, this class is responsible for differentiation of functions 
* This class is an abstract class defining an abstract method jacobian
* the jacobian method return the jacobian of a given vector function
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* Special Cases (ordered from general to specefic):
* + n = m
* In this case the method returns a square matrix, and there is a jabobian_det function which returns the determinant of the jacobian
* 
* + n = 1
* In this case the function is a vector univariate function:
* The method returns a column matrix, there is a derivative method which convert the column matrix as a vector  
* 
* + m = 1
* In this case the function is a scalar multivariate function:
* The method returns a row matrix, there is a gradient method which converts the column matrix as a vector
* 
* + n = m = 1
* In this case the function is a scalar univariate function
* The method returns a 1x1 matrix, there is a derivative method which converts the 1x1 matrix to a scalar
*/

namespace math_rz {

	template<typename F, int m, int n, typename E1, typename E2>
	class derivator
	{
	public:
		derivator(E1 p0) :x0(p0) {}
		virtual matrix<F, n, m> jacobian(const function<E1, E2>& f) const = 0;
	protected:
		E1 x0;
	};


	template<typename F, int n, typename E1, typename E2>
	class derivator<F,n,n,E1,E2>
	{
	public:
		derivator(E1 p0) :x0(p0) {}

		virtual square_matrix<F, n> jacobian(const function<E1, E2>& f) const = 0;
		virtual F jacobian_det(const function<E1, E2>& f) const
		{
			return jacobian(f).det();
		}
	protected:
		E1 x0;
	};


	template<typename F, int n, typename E2>
	class derivator<F, 1, n, F, E2>
	{
	public:
		derivator(F p0) :x0(p0) {}

		virtual matrix<F,n,1> jacobian(const function<F, E2>& f) const = 0;
		virtual E2 derivative(const function<F, E2>& f) const
		{
			E2 D;
			auto J = jacobian(f);
			for (int i = 0; i < n; i++)
				D[i] = J[i][0];
			return D;
		}
	protected:
		F x0;
	};


	template<typename F, int m, typename E1>
	class derivator<F, m,1, E1, F>
	{
	public:
		derivator(E1 p0) :x0(p0) {}

		virtual matrix<F, 1, m> jacobian(const function<E1, F>& f) const = 0;
		virtual E1 gradient(const function<E1, F>& f) const
		{
			E1 D;
			auto J = jacobian(f).transpose();
			for (int i = 0; i < m; i++)
				D[i] = J[i][0];
			return D;
		}
	protected:
		E1 x0;
	};

	template<typename F>
	class derivator<F, 1, 1, F, F>
	{
	public:
		derivator(F p0) :x0(p0) {}

		virtual square_matrix<F, 1> jacobian(const function<F, F>& f) const = 0;
		virtual F jacobian_det(const function<F, F>& f) const
		{
			return jacobian(f)[0][0];
		}
		virtual F gradient(const function<F, F>& f) const
		{
			return jacobian(f)[0][0];
		}
		virtual F derivative(const function<F, F>& f) const
		{
			return jacobian(f)[0][0];
		}
	protected:
		F x0;
	};

}