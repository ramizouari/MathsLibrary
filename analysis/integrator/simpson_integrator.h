#pragma once
#pragma once
#include "integrator.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"

namespace math_rz
{
	template<typename E, typename F>
	class simpsons_integrator :public  integrator<E, F>
	{
		integer cuts;
		real_field a, b;
	public:
		simpsons_integrator(real_field _a, real_field _b, integer _cuts = 100) :a(_a), b(_b), cuts(_cuts)
		{

		}
		F integrate(const function<E, F>& f) const override
		{
			F R;
			real_field u = std::min(a, b), v = std::max(a, b);
			real_field eps = (v - u) / cuts;
			int i = 0;
			for (real_field k1 = u,k2=k1+eps; k2 <= v; k1 += eps,k2+=eps,i++)
					R += (eps/2. ) * (f(k1)+2.*f((k1+k2)/2.)+f(k2));
			if (a < b)
				return R;
			else return -R;
		}
	};
}