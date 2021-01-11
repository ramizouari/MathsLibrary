#pragma once
#include "metric.h"
namespace math_rz::linalg::structure::vector
{
	template<typename F,int n>
	class norm_topology :public metric_topology<F,n>
	{
	protected:
		using vector_type = metric_topology<F, n>::vector_type;
	public:
		virtual real_field norm(const vector_type& a) const = 0;
		real_field metric(const vector_type& p, const vector_type& q) const
		{
			return  norm(p - q);
		}
	};

/*	template<typename F, int n, int m>
	class Linf_function_norm :public norm_topology<F, n,  m>
	{
		real_field a, b;
		integer cuts;
		
	public:
		Linf_norm(real_field _a, real_field _b, integer _c) :a(_a), b(_b), cuts(_c) {}
		real_field norm(const polynomial<F>& f) const
		{
			real_field eps = (b - a) / cuts;
			real_field u = a, R = 0;
			for (int i = 0; i <= cuts; i++, u += eps)
				R = std::max({ R,f(u).abs() });
			return R;
		}
	};*/
/*
	template<typename F, int n, int m>
	class L1_function_norm :public norm_topology<F, n, m>
	{
		std::shared_ptr<integrator<F, F>> I_ptr;
	public:
		L1_function_norm(std::shared_ptr<integrator<F, F>> _I_ptr) :I_ptr(_I_ptr) {}
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
	};*/
	/*
	template<typename F, int n, int m>
	class Lp_function_norm :public norm_topology<F, n, m>
	{
		integer p;
		std::shared_ptr<integrator<F, F>> I_ptr;
	public:
		Lp_function_norm(integer _p, std::shared_ptr<integrator<F, F>> _I_ptr) :I_ptr(_I_ptr), p(_p) {}
		Lp_function_norm(integer _p, integrator<F, F>* _I_ptr) :I_ptr(_I_ptr), p(_p) {}

		real_field norm(const polynomial<F>& f) const
		{
			return std::pow(I_ptr->integrate
			(
				general_function<F, F>([&](const F& u)->F
					{
						return std::pow(f(u).abs(), p);
					})
			), 1. / p);
		}
	};*/

	template<typename F, int n>
	class Lp_vector_norm :public norm_topology<F,n>
	{
		using vector_type = norm_topology<F, n>::vector_type;
		real_field p;
	public:
		Lp_vector_norm(real_field _p) :p(_p) {}

		real_field norm(const vector_type& M) const
		{
			real_field R = 0;
			constexpr auto inf = std::numeric_limits<long double>::infinity();
			for (int i = 0; i < n; i++)
				if (p < inf)
					R += std::pow(M[i].abs(), p);
				else R = std::max(R, M[i]);
			if (p < inf)
				R = std::pow(R, 1. / p);
			return R;
		}
	};
}