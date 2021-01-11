#pragma once
#include "integrator.h"
#include "analysis/function.h"
#include "analysis/normed_finite_dimensional_space.h"
#include <deque>

namespace math_rz {
	/*
	* This class is for integration of high dimensional (n>=2) Lp spaces
	* In pratice, due to exponential complexity, we advice that n<=10
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
	

		/*
		* This function calculates the multiple integral by the multi dimensional "rectangular method"
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



	/*
	* This class is for integration of high dimensional (n>=2) Lp spaces
	* In pratice, due to exponential complexity, we advice that n<=10
	*/
/*	template<typename F, int n, int p = 2>
	class stieltjes_multiple_integrator :public integrator <Lp_finite_dimensional_space< real_field, p, n>, F>
	{
		using Lp = Lp_finite_dimensional_space< real_field, p, n>;
		std::vector<integer> cuts_set;
		std::vector<real_field> lower_bounds, upper_bounds;
		function<Lp, F>& g;
	public:
		stieltjes_multiple_integrator(function<Lp,F>& _g,const std::vector<real_field>& L, const std::vector<real_field>& U, const std::vector<integer>& C)
			:lower_bounds(L), upper_bounds(U), cuts_set(C),g(_g)
		{}

		F integrate(const function<Lp, F>& f)const override
		{
			Lp w(lower_bounds);
			Lp eps;
			real_field V = 1;
			for (int i = 0; i < n; i++)
			{
				eps[i] = (upper_bounds[i] - lower_bounds[i]) / cuts_set[i];
			}
			std::vector<integer> counts(n, 0);
			counts.reserve(n);
			F R;
			int m = counts.size();
			while (m > 0)
			{
				for (int i = m; i < n; i++)
				{
					counts[i] = 0;
					w[i] = lower_bounds[i];
				}
				m = n;
				//using result_type = std::conditional_t<F::dimension == 1, std::vector<std::vector<real_field>>, real_field>;
				std::vector<std::vector<real_field>> M;
				for (int i = 0; i < n; i++)
				{
					F g1 = g(w);
					w[i] += eps[i];
					F g2 = g(w);
					if constexpr (F::dimension > 1)
						M.push_back((g2 - g1).get_vect());
					else M.push_back({ g2 - g1 });
					w[i] -= eps[i];
				}
				square_matrix<real_field, F::dimension> W(std::move(M));
				R += W.det() * f(w);
				while (m > 0 && counts[m - 1] == (cuts_set[m - 1] - 1))
				{
					m--;
				}
				if (m > 0)
				{
					counts[m - 1]++;
					w[m - 1] += eps[m - 1];
				}
			}
			return R;
		}
	};*/
}