#pragma once
#include "function.h"
#include "real_field.h"
#include "integrator.h"

class rr_integrable_function :public function<real_field,real_field>
{
public:
	real_field integral(const integrator<real_field,real_field> &w) const
	{
		return w.integrate(*this);
	}
};