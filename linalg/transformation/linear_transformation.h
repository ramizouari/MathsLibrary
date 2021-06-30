#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"

namespace math_rz::linalg
{
	/*
	* This class is reserved for linear transformations which can be calculated in o(n²)
	*/
	template<vector_space_constraint::vector_space E,vector_space_constraint::vector_space F>
	class linear_transformation:public analysis::function<E,F>
	{
		template<int m>
		using matrix_type = matrix<typename E::base_field, E::dimension, m>;
	public:
		template<int m>
		matrix_type<m> operator()(const matrix_type<m>& A) const
		{
			std::vector<std::vector<typename E::base_field>> B;
			for (int i = 0; i < m; i++)
				B.push_back(this->operator()(A.get_column(i)).get_vect());
			return B;
		}
	};


	template<vector_space_constraint::vector_space E>
	class endomorphism :public linear_transformation<E, E>
	{
		template<int m>
		using matrix_type = matrix<typename E::base_field, E::dimension,m>;
	public:
		virtual E& apply(E& u) const = 0;
		template<int m>
		matrix_type<m>& apply(matrix_type<m>& A) const
		{
			return A = std::move(this->operator()(A));
		}
		
		E operator()(const E& u) const
		{
			E v(u);
			return apply(v);
		}
	};
}