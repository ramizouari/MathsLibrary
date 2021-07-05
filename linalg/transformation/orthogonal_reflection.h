#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linear_transformation.h"

namespace math_rz::linalg
{
	template<vector_space_constraint::inner_product_space E>
	class orthogonal_reflection :public endomorphism<E>
	{
		using K = typename E::base_field;
		inline static constexpr int n = E::dimension;
		std::vector<E> V;
	public:
		orthogonal_reflection(std::vector<E> _V) :V(_V)
		{
			for (auto& v : V)
				v /= K(v.norm());
		}
		orthogonal_reflection(E s) :orthogonal_reflection(std::vector<E>({ s })) {}
		E& apply(E& u) const
		{
			E R;
			for (const auto& v : V)
			{
				K A = K(2) * v.inner_product(u);
				R += A* v;
			}
			return u -=R;
		}

		bool is_zero() const
		{
			return false;
		}
	};
}