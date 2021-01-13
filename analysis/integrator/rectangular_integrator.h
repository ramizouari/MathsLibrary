#pragma once
#include "integrator.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"

namespace math_rz::analysis
{
	template<typename E,typename F>
	class rectangular_integrator :public  integrator<E, F>
	{
		integer cuts;
		real_field a, b;
	public:
		rectangular_integrator(real_field _a, real_field _b, integer _cuts = 100) :a(_a), b(_b), cuts(_cuts)
		{

		}
		F integrate(const function<E, F>& f) const override
		{
			real_field eps = (b - a) / cuts;
			real_field u = a;
			F result;
			for (int i = 0; i < cuts; i++, u += eps)
				result += eps * f(u);
			return result;
		}
	};


	/*template<int n,int m >
	class default_integrator :public  integrator<coordinate_space<real_field, n>,coordinate_space<real_field,m> >
	{
		using E1 = coordinate_space<real_field, n>;
		using E2 = coordinate_space<real_field, m>;
		integer cuts;
		E1 a, b;
	public:
		default_integrator(real_field _a, real_field _b, integer _cuts = 100) :a(_a), b(_b), cuts(_cuts)
		{

		}
		F integrate(const function<E1, E2>& f) const override
		{
			F R;
			real_field u = std::min(a, b), v = std::max(a, b);
			real_field eps = (v - u) / cuts;
			for (real_field k = u; k <= v; k += eps)
				R += eps * f(k);
			if (a < b)
				return R;
			else return -R;
		}
	};*/
}