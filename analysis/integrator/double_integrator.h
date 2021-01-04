#pragma once
#include "integrator.h"
#include "analysis/normed_finite_dimensional_space.h"


namespace math_rz {
	template<typename F,int p=2>
	class double_integrator :public integrator <Lp_finite_dimensional_space< real_field,p,2>,F>
	{
		using Lp = Lp_finite_dimensional_space< real_field, p, 2>;
		integer cuts1, cuts2;
		real_field a, b, c, d;
	public:
		double_integrator(real_field _a,real_field _b,real_field _c,real_field _d,integer _cuts1,integer _cuts2)
			:a(_a),b(_b),c(_c),d(_d),cuts1(_cuts1),cuts2(_cuts2)
		{}

		F integrate(const function<Lp, F>& f)const override
		{
			real_field eps2 = (d - c) / cuts2, eps1 = (b - a) / cuts1;
			F result;
			real_field u = a, v = c;
			for (int i = 0; i < cuts1; i++, v = c, u += eps1)
				for (int j = 0; j < cuts2; j++, v += eps2)
					result += (eps1 * eps2) *f(Lp({ u,v }));
			return result;
		}
	};


	template<typename F, int p = 2>
	class trapezoidal_double_integrator :public integrator <Lp_finite_dimensional_space< real_field, p, 2>, F>
	{
		using Lp = Lp_finite_dimensional_space< real_field, p, 2>;
		integer cuts1, cuts2;
		real_field a, b, c, d;
	public:
		trapezoidal_double_integrator(real_field _a, real_field _b, real_field _c, real_field _d, integer _cuts1, integer _cuts2)
			:a(_a), b(_b), c(_c), d(_d), cuts1(_cuts1), cuts2(_cuts2)
		{}

		F integrate(const function<Lp, F>& f) const override
		{
			real_field eps2 = (d - c) / cuts2, eps1 = (b - a) / cuts1;
			F result;
			real_field u=a, v=c;
			for (int i=0;i<cuts1;i++,v=c,u+=eps1)
				for (int j=0;j<cuts2;j++,v+=eps2)
					result += (eps1 * eps2*real_field(.25)) * 
					(
						f(Lp({ u,v })) +
						f(Lp({ u + eps1,v })) +
						f(Lp({ u,v+eps2 })) +
						f(Lp({ u+eps1,v+eps2 }))
					);
			return result;
		}
	};
}