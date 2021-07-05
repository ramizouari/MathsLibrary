#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linear_transformation.h"
#include "orthogonal_reflection.h"

namespace math_rz::linalg
{
	template<vector_space_constraint::inner_product_space E>
	class house_holder_reflection :public  orthogonal_reflection<E>
	{
		using K = typename E::base_field;
	public:
		house_holder_reflection(E s) :orthogonal_reflection<E>({ s / K(s.norm()) }) {}
		house_holder_reflection(std::vector<K> s) :house_holder_reflection(E(s)) {}
	};
}