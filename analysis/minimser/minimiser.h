#pragma once
#include "analysis/function.h"
#include "real_field.h"
namespace math_rz::analysis
{
	template<typename E>
	class minimiser
	{
	public:
		virtual E argmin(const function<E, real_field>& f) const = 0;
		virtual real_field min(const function<E, real_field>& f) const
		{
			return f(argmin(f));
		}
	};
}