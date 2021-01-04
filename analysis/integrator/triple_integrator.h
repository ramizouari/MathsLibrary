#pragma once
#include "integrator.h"
#include "analysis/normed_finite_dimensional_space.h"


namespace math_rz {
	template<typename F, int p = 2>
	class triple_integrator :public integrator <Lp_finite_dimensional_space< real_field, p, 3>, F>
	{
		using Lp = Lp_finite_dimensional_space< real_field, p, 3>;
		integer cuts1, cuts2,cuts3;
		real_field a, b, c, d,e,f;
	public:
		triple_integrator(real_field _a, real_field _b, real_field _c, real_field _d, 
			real_field _e,real_field _f, integer _cuts1, integer _cuts2,integer _cuts3)
			:a(_a), b(_b), c(_c), d(_d),e(_e),f(_f), cuts1(_cuts1), cuts2(_cuts2),cuts3(_cuts3)
		{}

		F integrate(const function<Lp, F>& h)const override
		{
			real_field eps3=(f-e)/cuts3,eps2 = (d - c) / cuts2, eps1 = (b - a) / cuts1;
			real_field result = 0;
			real_field u, v, w;
			for (int i = 0; i < cuts1; i++, v = c, u += eps1,v=0,w=0)
				for (int j = 0; j < cuts2; j++, v += eps2,w=0)
					for(int k=0;k<cuts3;k++,w+=eps3)
					result += (eps1 * eps2*eps3) * h(Lp({ u,v,w }));
			return result;
		}
	};


	template<typename F, int p = 2>
	class trapezoidal_triple_integrator :public integrator <Lp_finite_dimensional_space< real_field, p, 3>, F>
	{
		using Lp = Lp_finite_dimensional_space< real_field, p, 3>;
		integer cuts1, cuts2, cuts3;
		real_field a, b, c, d, e, f;
	public:
		trapezoidal_triple_integrator(real_field _a, real_field _b, real_field _c, real_field _d,
			real_field _e, real_field _f, integer _cuts1, integer _cuts2, integer _cuts3)
			:a(_a), b(_b), c(_c), d(_d), e(_e), f(_f), cuts1(_cuts1), cuts2(_cuts2), cuts3(_cuts3)
		{}

		F integrate(const function<Lp, F>& h)const override
		{
			real_field eps3 = (f - e) / cuts3, eps2 = (d - c) / cuts2, eps1 = (b - a) / cuts1;
			real_field result = 0;
			real_field u, v, w;
			for (int i = 0; i < cuts1; i++, v = c, u += eps1, v = 0, w = 0)
				for (int j = 0; j < cuts2; j++, v += eps2, w = 0)
					for (int k = 0; k < cuts3; k++, w += eps3)
						result += (eps1 * eps2 * eps3*real_field(.125)) * 
						(
							h(Lp({ u,v,w }))+ h(Lp({ u+eps1,v,w }))+
							h(Lp({ u,v+eps2,w })) + h(Lp({ u + eps1,v+eps2,w })) +
							h(Lp({ u,v+eps2,w+eps3 })) + h(Lp({ u + eps1,v+eps2,w+eps3 })) +
							h(Lp({ u,v,w+eps3 })) + h(Lp({ u + eps1,v,w+eps3 }))
						);
			return result;
		}
	};
}