#pragma once
#include "integrator.h"
#include "analysis/normed_finite_dimensional_space.h"
#include <deque>

namespace math_rz {
	/*
	* This class is for integration of high dimensional (n>=2) Lp spaces
	* In pratice,because of complexity, we advice that n<=10
	*/
	template<typename F,int n, int p = 2>
	class multiple_integrator :public integrator <Lp_finite_dimensional_space< real_field, p, n>, F>
	{
		using Lp = Lp_finite_dimensional_space< real_field, p, n>;
		std::vector<integer> cuts_set;
		std::vector<real_field> lower_bounds, upper_bounds;
	public:
		multiple_integrator(const std::vector<real_field>& L, const std::vector<real_field>& U, const std::vector<integer> &C)
			:lower_bounds(L),upper_bounds(U),cuts_set(C)
		{}
		multiple_integrator(const real_field& L, const std::vector<real_field>& U, const std::vector<integer>& C)
			:lower_bounds(n,L), upper_bounds(U), cuts_set(C)
		{}

		multiple_integrator(const std::vector<real_field>& L, const real_field& U, const std::vector<integer>& C)
			:lower_bounds(L), upper_bounds(n,U), cuts_set(C)
		{}
		multiple_integrator(const real_field& L, const real_field& U, const std::vector<integer>& C)
			:lower_bounds(n,L), upper_bounds(n,U), cuts_set(C)
		{}


		multiple_integrator(const std::vector<real_field>& L, const std::vector<real_field>& U, 
			const integer& C)
			:lower_bounds(L), upper_bounds(U), cuts_set(n,C)
		{}
		multiple_integrator(const real_field& L, const std::vector<real_field>& U, const integer& C)
			:lower_bounds(n,L), upper_bounds(U), cuts_set(n,C)
		{}

		multiple_integrator(const std::vector<real_field>& L, const real_field& U, const integer& C)
			:lower_bounds(L), upper_bounds(n,U), cuts_set(n,C)
		{}
		multiple_integrator(const real_field& L, const real_field& U, const integer& C)
			:lower_bounds(n,L), upper_bounds(n,U), cuts_set(n,C)
		{}


		/*
		* This function calculates the multiple integral by the multi dimensional "rectangle method"
		* Let W the product of elements of cuts_set
		* The complexity is O(W)
		*/
		F integrate(const function<Lp, F>& f)const override
		{
			Lp w(lower_bounds);
			Lp eps;
			real_field V=1;
			for (int i = 0; i < n; i++)
			{
				eps[i] = (upper_bounds[i] - lower_bounds[i]) / cuts_set[i];
				V *= eps[i];
			}
			std::vector<integer> counts(n,0);
			counts.reserve(n);
			F R;
			int m = counts.size();
			while (m>0)
			{
				for (int i = m; i < n; i++)
				{
					counts[i]=0;
					w[i] = lower_bounds[i];
				}
				m = n;
				R += V*f(w);
				while (m>0 && counts[m-1] == (cuts_set[m - 1]-1))
				{
					m--;
				}
				if (m>0)
				{
					counts[m-1]++;
					w[m - 1] += eps[m - 1];
				}
			}
			return R;
		}
	};
}