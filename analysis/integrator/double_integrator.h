#pragma once
#include "integrator.h"
#include "analysis/normed_finite_dimensional_space.h"
/*
* This is a 
*/

namespace math_rz::analysis {
	template <linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F > requires (E::dimension == 2)
	class double_integrator :public integrator <E,F>
	{
		using K = typename E::base_field;
		integer cuts1, cuts2;
		K a, b, c, d;
	public:
		double_integrator(K _a,K _b,K _c,K _d,integer _cuts1,integer _cuts2)
			:a(_a),b(_b),c(_c),d(_d),cuts1(_cuts1),cuts2(_cuts2)
		{}

		F integrate(const function<E, F>& f)const override
		{
			K eps2 = (d - c) / cuts2, eps1 = (b - a) / cuts1;
			F result;
			K u = a, v = c;
			for (int i = 0; i < cuts1; i++, v = c, u += eps1)
				for (int j = 0; j < cuts2; j++, v += eps2)
				{
					E s({ u,v });
					result += (eps1 * eps2) * f(s);
				}
			return result;
		}
	};


	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F > requires (E::dimension == 2)
	class trapezoidal_double_integrator :public integrator <E, F>
	{
		using K = typename E::base_field;
		integer cuts1, cuts2;
		K a, b, c, d;
	public:
		trapezoidal_double_integrator(K _a, K _b, K _c, K _d, integer _cuts1, integer _cuts2)
			:a(_a), b(_b), c(_c), d(_d), cuts1(_cuts1), cuts2(_cuts2)
		{}

		F integrate(const function<E, F>& f) const override
		{
			K eps2 = (d - c) / K(cuts2), eps1 = (b - a) / K(cuts1);
			F result;
			K u=a, v=c;
			for (int i=0;i<cuts1;i++,v=c,u+=eps1)
				for (int j=0;j<cuts2;j++,v+=eps2)
					result += K(eps1 * eps2*real_field(.25)) *
					(
						f(E({ u,v })) +
						f(E({ u + eps1,v })) +
						f(E({ u,v+eps2 })) +
						f(E({ u+eps1,v+eps2 }))
					);
			return result;
		}
	};

	template<linalg::vector_space_constraint::vector_space E, linalg::vector_space_constraint::vector_space F > requires (E::dimension == 2)
		class simpson_double_integrator :public integrator <E, F>
	{
		using K = typename E::base_field;
		integer cuts1, cuts2;
		K a, b, c, d;
	public:
		simpson_double_integrator(K _a, K _b, K _c, K _d, integer _cuts1, integer _cuts2)
			:a(_a), b(_b), c(_c), d(_d), cuts1(_cuts1), cuts2(_cuts2)
		{}

		F integrate(const function<E, F>& f) const override
		{
			K eps2 = (d - c) / K(cuts2), eps1 = (b - a) / K(cuts1);
			F result;
			K u = a, v = c;
			for (int i = 0; i < cuts1; i++, v = c, u += eps1)
				for (int j = 0; j < cuts2; j++, v += eps2)
					result += K(eps1 * eps2 /36) *
					(
						f(E({ u,v })) +
						4*f(E({ u + eps1/2,v })) +
						f(E({ u + eps1,v }))+
						4*f(E({ u,v + eps2/2 })) +
						16*f(E({ u + eps1/2,v + eps2/2 }))+
						4*f(E({ u,v+eps2/2 })) +
						f(E({ u+ eps1 ,v })) +
						4 * f(E({ u + eps1,v + eps2 / 2 })) +
						f(E({ u+ eps1,v + eps2 })) 
						);
			return result;
		}
	};
}