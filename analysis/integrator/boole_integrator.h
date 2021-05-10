#pragma once
#pragma once
#include "integrator.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"

namespace math_rz::analysis
{
	template<typename E, typename F>
	class boole_integrator :public  integrator<E, F>
	{
		integer cuts;
		real_field a, b;
	public:
		boole_integrator(real_field _a, real_field _b, integer _cuts = 100) :a(_a), b(_b), cuts(_cuts)
		{

		}
		F integrate(const function<E, F>& f) const override
		{
			F R = 0;
			real_field u = std::min(a, b), v = std::max(a, b);
			real_field eps = (v - u) / cuts;
			int i = 0;
			for (real_field k1 = u, k2 = k1 + eps; k2 <= v; k1 += eps, k2 += eps, i++)
				R += (eps / 90) * (7.*f(k1) + 32. * f((3.*k1 + k2) / 4.) + 
					12. * f((k1 + k2)/2.) + 32. * f((k1 + 3.*k2)/4.) + 7.*f(k2));
			if (a < b)
				return R;
			else return -R;
		}
	};
}