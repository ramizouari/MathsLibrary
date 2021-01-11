#pragma once
#include "analysis/function.h"
#include "real_field.h"
#include "integer.h"

/*
* This is the base class of integrators
* - Informally speaking, an integrator is an object which gives "an integral" of a function from a set A to a set B
* - An integrator is an operator from the vector space of functions to another normed vector space 
*/

namespace math_rz {
	template <typename A, typename B,typename I=B>
	class integrator
	{
	public:
		using domain = A;
		using codomain = B;
		using result_type = I;
		virtual I integrate(const function<A, B>& f) const = 0;
	};

	/*
	* This class gives an implementation of trapezoidal rule of a univariate (possibly vector) function
	*/

	template<typename E,typename F>
	class trapezoidal_integrator :public  integrator<E, F,F>
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
			real_field u = a;
			for (int i = 0; i < cuts; i++, u += eps)
				R += real_field(eps._v*.5) * (f(u)+f(u+eps));
			return R;
		}
	};
}