#pragma once
#include "minimiser.h"
namespace math_rz::analysis
{
	template<typename E>
	class constrained_minimiser:public minimiser<E>
	{
	public:
		virtual E argmin(const function<E, real_field>& f) const = 0;
		virtual real_field min(const function<E, real_field>& f) const
		{
			return f(argmin(f));
		}
	};
}