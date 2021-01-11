#pragma once
#include "real_field.h"
#include "poly/polynomial.h"
#include "metric.h"
#include "analysis/integrator/integrator.h"
namespace math_rz::poly::structure
{
	template<typename F>
	class norm_topology :public metric_topology<F>
	{
	public:
		virtual real_field norm(const polynomial<F>& a) const = 0;
		real_field metric(const free_algebra<F>& p, const free_algebra<F>& q) const
		{
			return  norm(p - q);
		}
	};

	template<typename F>
	class Linf_function_norm :public norm_topology<F>
	{
		real_field a, b;
		integer cuts;
	public:
		Linf_function_norm(real_field _a,real_field _b,integer _c):a(_a),b(_b),cuts(_c){}
		real_field norm(const polynomial<F>& f) const
		{
			real_field eps = (b - a) / cuts;
			real_field u = a,R=0;
			for (int i = 0; i <= cuts; i++, u += eps)
				R = std::max({ R,f(u).abs() });
			return R;
		}
	};

	template<typename F>
	class L1_function_norm :public norm_topology<F>
	{
		std::shared_ptr<integrator<F, F>> I_ptr;
	public:
		L1_function_norm(std::shared_ptr<integrator<F,F>> _I_ptr) :I_ptr(_I_ptr) {}
		L1_function_norm(integrator<F, F>* _I_ptr) :I_ptr(_I_ptr) {}
		real_field norm(const polynomial<F>& f) const
		{
			return I_ptr->integrate
			(
				general_function<F, F>([&](const F& u)->F
					{
						return f(u).abs();
					})
			);
		}
	};

	template<typename F>
	class Lp_function_norm :public norm_topology<F>
	{
		real_field p;
		std::shared_ptr<integrator<F, F>> I_ptr;
	public:
		Lp_function_norm(real_field _p,std::shared_ptr<integrator<F, F>> _I_ptr) :I_ptr(_I_ptr), p(_p) {}
		Lp_function_norm(real_field _p,integrator<F, F>* _I_ptr) :I_ptr(_I_ptr), p(_p) {}

		real_field norm(const polynomial<F>& f) const
		{
			constexpr auto inf = std::numeric_limits<long double>::infinity();
			if (p<inf)
				return std::pow(I_ptr->integrate
				(
					general_function<F, F>([&](const F& u)->F
						{
							return std::pow(f(u).abs(),p);
						})
				),1./p);
			else
			{
				// To Do: change this
				return Linf_function_norm<F>(-1, 1, 100).norm();
			}
		}
	};

	template<typename F>
	class Lp_vector_norm :public norm_topology<F>
	{
		real_field p;
	public:
		Lp_vector_norm(real_field _p):p(_p) {}

		real_field norm(const polynomial<F>& f) const
		{
			constexpr auto inf = std::numeric_limits<long double>::infinity();
			real_field R = 0;
			for (const auto& v : f.get_vect())
				if (p < inf)
					R += std::pow(v.abs(), p);
				else R = std::max(R, v);
			if (p < inf)
				R = std::pow(R, 1. / p);
			return R;
		}
	};
}