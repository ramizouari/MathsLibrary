#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linear_transformation.h"

namespace math_rz::linalg
{
	template<vector_space_constraint::inner_product_space E>
	class orthogonal_projection :public endomorphism<E>
	{
		using K = typename E::base_field;
		inline static constexpr int n = E::dimension;
		std::vector<E> V;
	public:
		orthogonal_projection(std::vector<E> _V) :V(_V) 
		{
			for (auto& v : V)
				v /= K(v.norm());
		}
		orthogonal_projection(E s) :orthogonal_projection(std::vector<E>({s})) {}
		E& apply(E& u) const
		{
			E R = u;
			for (const auto &v:V)
				R += v.inner_product(u);
			return u = std::move(R);
		}

		bool is_zero() const
		{
			return false;
		}
	};
}