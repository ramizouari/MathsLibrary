#pragma once
#include "analysis/function.h"
#include "real_field.h"
namespace math_rz::analysis
{
	template<typename E>
	class root_finder
	{
	public:
		virtual E root(const function<E, real_field>& f) const = 0;
	};
}