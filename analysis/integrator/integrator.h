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
	class trapezoidal_integrator :public  integrator<E, F>
	{
		integer cuts;
		real_field a, b;
	public:
		trapezoidal_integrator(real_field _a, real_field _b, integer _cuts = 100) :a(_a), b(_b), cuts(_cuts)
		{

		}
		F integrate(const function<E, F>& f) const override
		{
			F R;
			real_field eps = (b - a) / cuts;
			real_field u = a,result;
			for (int i = 0; i < cuts; i++, u += eps)
				result += (eps*real_field(.5)) * (f(u)+f(u+eps));
			return result;
		}
	};
}