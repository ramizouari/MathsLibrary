#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linear_transformation.h"

namespace math_rz::linalg
{
	template<vector_space_constraint::inner_product_space E>
	class rotation :public endomorphism<E>
	{
		using K = typename E::base_field;
		inline static constexpr int n = E::dimension;
		E u,v;
		K a, b;
	public:
		rotation(E _u,E _v,K _a,K _b) :u(_u),v(_v),a(_a),b(_b)
		{
			v /= K(v.norm());
			u /= K(u.norm());
			K r = std::sqrt(pow(a.norm(), 2) + pow(b.norm(), 2));
			a /= r;
			b /= r;
		}
		E& apply(E& w) const
		{
			auto s=u.inner_product(w), t=v.inner_product(w);
			K p = (a * s - b * t - K(1)),q=(b*s+a*t-K(1));
			E R = p*u+q*v;
			return w += R;
		}

		bool is_zero() const
		{
			return false;
		}
	};
}