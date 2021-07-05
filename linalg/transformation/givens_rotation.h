#pragma once
#include "analysis/function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linear_transformation.h"

namespace math_rz::linalg
{
	template<vector_space_constraint::inner_product_space E>
	class givens_rotation :public endomorphism<E>
	{
		using K = typename E::base_field;
		int p, q;
		K a, b;
	public:
		givens_rotation(int _p,int _q,K _a,K _b) :p(_p),q(_q) 
		{
			K r = std::sqrt(pow(_a.norm(), 2) + pow(_b.norm(), 2));
			a = _a / r;
			b = _b / r;
		}

		givens_rotation(int _p, int _q, K theta) :p(_p), q(_q),a(std::cos(theta)),b(std::sin(theta))
		{
		}
		E& apply(E& u) const
		{
			K s = u[p], t = u[q];
			u[p] = a * s - b * t;
			u[q] = b * s + a * t;
			return u;
		}

		bool is_zero() const
		{
			return false;
		}
	};
}