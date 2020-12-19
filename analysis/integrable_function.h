#pragma once
#include "function.h"
#include "real_field.h"
#include "integrator.h"

namespace math_rz {
	template<typename F>
	class curve :public function<real_field, F>
	{
	public:
		F integral(const integrator<real_field, F>& w) const
		{
			return w.integrate(*this);
		}
	};

	template<typename F>
	class complex_curve :public function<complex, F>
	{

	};
}