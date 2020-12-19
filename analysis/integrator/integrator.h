#pragma once
#include "analysis/function.h"
#include "real_field.h"
#include "integer.h"

namespace math_rz {
	template <typename A, typename B>
	class integrator
	{
	public:
		virtual B integrate(const function<A, B>& f) const = 0;
	};



	template<typename E,typename F>
	class trapezoid_integrator :public  integrator<E, F>
	{
		integer cuts;
		real_field a, b;
	public:
		trapezoid_integrator(real_field _a, real_field _b, integer _cuts = 100) :a(_a), b(_b), cuts(_cuts)
		{

		}
		F integrate(const function<E, F>& f) const override
		{
			F R;
			real_field u = std::min(a, b), v = std::max(a, b);
			real_field eps = (v - u) / cuts;
			for (real_field k = u; k < v; k += eps)
			{
				real_field s = eps / 2;
				R += s * (f(k) + f(k + eps));
			}
			if (a < b)
				return R;
			else return -R;
		}
	};
}