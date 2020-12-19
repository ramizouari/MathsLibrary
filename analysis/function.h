#pragma once
#include "linalg/vector_space.h"

namespace math_rz {
	template<typename A, typename B>
	class abstract_function
	{
	public:
		virtual B operator()(const A& a) const = 0;

	};


	template<typename A, typename B>
	class summable_function : public abstract_function<A, B>, public vector_space<B>
	{
	public:

	};

	template<typename A, typename B>
	class function : public summable_function<A, B>
	{
	};
}