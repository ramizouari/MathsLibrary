#pragma once
#include "linalg/vector_space.h"
#include "real_field.h"

namespace math_rz {
	template<typename F>
	class normed_space : public vector_space<F>
	{
		virtual real_field norm() const = 0;
	};
}