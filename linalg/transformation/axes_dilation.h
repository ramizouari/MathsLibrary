#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linear_transformation.h"

namespace math_rz::linalg
{
	template<vector_space_constraint::vector_space E>
	class axes_dilation :public endomorphism<E>
	{
		using K = typename E::base_field;
		inline static constexpr int n = E::dimension;
		K k;
		E lambda;
	public:
		axes_dilation(E s) :lambda(s) {}
		axes_dilation(std::vector<K> s) :axes_dilation(E(s)) {}
		E& apply(E& u) const
		{
			for (int i = 0; i < n; i++)
				u[i] *= lambda[i];
			return u;
		}

		bool is_zero() const
		{
			return lambda.is_zero();
		}
	};
}