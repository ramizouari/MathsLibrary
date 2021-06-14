#pragma once
#include "abstract_module.h"

namespace math_rz {
	template<typename F>
	class vector_space :virtual public abstract_module<F>
	{
	public:
		using base_field = F;
	};
}

