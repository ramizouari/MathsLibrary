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

	};


	template<vector_space_constraint::vector_space E>
	class endomorphism :public linear_transformation<E, E>
	{
	public:
		virtual E& apply(E& u) const = 0;
		E operator()(const E& u) const
		{
			E v(u);
			return apply(v);
		}
	};
}