#pragma once
#include "function.h"
#include "real_field.h"
#include "integer.h"

template <typename A,typename B>
class integrator
{
public:
	virtual B integrate(const function<A, B>& f) const = 0;
};

class rectangle_integrator:public  integrator<real_field,real_field>
{
	integer cuts;
	real_field a, b;
public:
	rectangle_integrator(real_field _a, real_field _b, integer _cuts = 100):a(_a),b(_b),cuts(_cuts)
	{

	}
	real_field integrate(const function<real_field, real_field>& f) const override
	{
		real_field R,u=std::min(a,b),v=std::max(a,b);
		real_field eps = (v - u) / cuts;
		for (real_field k = u; k <= v; k += eps)
			R += f(k) * eps;
		if (a < b)
			return R;
		else return -R;
	}
};