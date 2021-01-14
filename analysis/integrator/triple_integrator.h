#pragma once
#include "integrator.h"
#include "analysis/normed_finite_dimensional_space.h"


namespace math_rz::analysis {
	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F> requires (E::dimension == 3)
	class triple_integrator :public integrator <E, F>
	{
		using K = typename E::base_field;
		integer cuts1, cuts2,cuts3;
		K a, b, c, d,e,f;
	public:
		triple_integrator(K _a, K _b, K _c, K _d, K _e,K _f, integer _cuts1, integer _cuts2,integer _cuts3)
			:a(_a), b(_b), c(_c), d(_d),e(_e),f(_f), cuts1(_cuts1), cuts2(_cuts2),cuts3(_cuts3)
		{}

		F integrate(const function<E, F>& h)const override
		{
			K eps3=(f-e)/K(cuts3),eps2 = (d - c) / K(cuts2), eps1 = (b - a) / K(cuts1);
			F result;
			//u,v,w are dummy variables used for integration
			K u=e, v=c, w=a;
			//This nested for loop of depth 3 calculates the integral
			for (int i = 0; i < cuts1; i++, v = c, u += eps1,v=c,w=a)
				for (int j = 0; j < cuts2; j++, v += eps2,w=a)
					for(int k=0;k<cuts3;k++,w+=eps3)
						result += (eps1 * eps2*eps3) * h(E({ u,v,w }));
			return result;
		}
	};


	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F> requires (E::dimension == 3)
	class trapezoidal_triple_integrator :public integrator <E, F>
	{
		using K = typename E::base_field;
		integer cuts1, cuts2, cuts3;
		K a, b, c, d, e, f;
	public:
		trapezoidal_triple_integrator(K _a, K _b, K _c, K _d,
			K _e, K _f, integer _cuts1, integer _cuts2, integer _cuts3)
			:a(_a), b(_b), c(_c), d(_d), e(_e), f(_f), cuts1(_cuts1), cuts2(_cuts2), cuts3(_cuts3)
		{}

		F integrate(const function<E, F>& h)const override
		{
			K eps3 = (f - e) / K(cuts3), eps2 = (d - c) / K(cuts2), eps1 = (b - a) / K(cuts1);
			F result;
			K u=e, v=c, w=a;
			for (int i = 0; i < cuts1; i++, v = c, u += eps1, v = c, w = a)
				for (int j = 0; j < cuts2; j++, v += eps2, w = a)
					for (int k = 0; k < cuts3; k++, w += eps3)
						result += (eps1 * eps2 * eps3*K(.125)) *
						(
							h(E({ u,v,w }))+ h(E({ u+eps1,v,w }))+
							h(E({ u,v+eps2,w })) + h(E({ u + eps1,v+eps2,w })) +
							h(E({ u,v+eps2,w+eps3 })) + h(E({ u + eps1,v+eps2,w+eps3 })) +
							h(E({ u,v,w+eps3 })) + h(E({ u + eps1,v,w+eps3 }))
						);
			return result;
		}
	};
}