#pragma once
#include "analysis/function.h"
#include "integrator.h"

namespace math_rz
{
	template<typename B>
	class integral_function : public function<real_field,B>
	{
		integral_function(function<real_field, B> base_function, integrator& _I) :f(base_function), I(_I)
		{}

		B operator()(const A& w) const
		{
			return 
		}
	private:
		function<real_field, B> f;
		integrator& I;
	};
}