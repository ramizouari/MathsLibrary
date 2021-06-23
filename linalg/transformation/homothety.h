#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linear_transformation.h"

namespace math_rz::linalg
{
	template<vector_space_constraint::vector_space E>
	class homothety :public endomorphism<E>
	{
		using K = typename E::base_field;
		K k;
	public:
		homothety(K a):k(a){}
		E& apply(E& u) const
		{
			u.foreach([k](auto& a) 
				{
					a *= k;
				});
			return u;
		}

		bool is_zero() const
		{
			return k.is_zero();
		}
	};
}