#pragma once
#include "linalg/matrix.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "analysis/function.h"
#include <type_traits>

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

namespace math_rz::analysis {

	template<typename E1, typename E2> 
	requires linalg::vector_space_constraint::vector_space_over_same_base_field<E1,E2>
	class derivator
	{
	public:
		inline static constexpr int domain_dimension = E1::dimension;
		inline static constexpr int codomain_dimension = E2::dimension;
		using base_field = typename E1::base_field;
		using matrix_type = std::conditional_t<domain_dimension == codomain_dimension,
			linalg::square_matrix<base_field, codomain_dimension>, 
			linalg::matrix<base_field,codomain_dimension, domain_dimension>>;

		virtual matrix_type jacobian(const function<E1, E2>& f,const E1& x0) const = 0;
		base_field jacobian_det(const function<E1, E2>& f, const E1& x0) const 
			requires (codomain_dimension==domain_dimension)
		{
			return jacobian(f, x0).det();
		}
		E2 derivative(const function<base_field, E2>& f, const base_field& x0) const
			requires (domain_dimension == 1)
		{
			return jacobian(f, x0).as_vector();
		}

		E1 gradient(const function<E1, base_field>& f, const E1& x0) const
			requires (codomain_dimension == 1)
		{
			return jacobian(f, x0).as_vector();
		}
	protected:
	};
}