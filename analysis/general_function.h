#pragma once
#include <functional>
#include "function.h"

namespace math_rz
{
	template<typename A,typename B>
	class general_function:public function<A,B>
	{
	public :
		general_function(std::function<B(A)> &&_f):f(std::move(_f)){}

		B operator()(const A& s) const
		{
			return f(s);
		}


		/*
		* This problem is undecidable, we will return just false with the unformal meaning:
		* almost all functions are not null
		*/
		bool is_zero() const override
		{
			return false;
		}
	private:
		std::function<B(A)> f;
	};
}