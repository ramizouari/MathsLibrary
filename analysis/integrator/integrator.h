#pragma once
#include "analysis/function.h"
#include "real_field.h"
#include "integer.h"

/*
* This is the base class of integrators
* - Informally speaking, an integrator is an object which gives "an integral" of a function from a set A to a set B
* - An integrator is an operator from the vector space of functions to another normed vector space 
*/

namespace math_rz::analysis {
	template <typename A, typename B,typename I>
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

	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F>
	class trapezoidal_integrator :public  integrator<E, F,F>
	{
		using K = typename E::base_field;
		integer cuts;
		K a, b;
	public:
		trapezoidal_integrator(K _a, K _b, integer _cuts = 100) :a(_a), b(_b), cuts(_cuts)
		{

		}
		F integrate(const function<E, F>& f) const override
		{
			F R;
			K eps = (b - a) / K(cuts);
			K u = a;
			for (int i = 0; i < cuts; i++, u += eps)
				R += (eps*K(.5)) * (f(u)+f(u+eps));
			return R;
		}
	};
}