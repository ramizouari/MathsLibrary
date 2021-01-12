#pragma once
#include "metric.h"
#include "linalg/finite_dimensional_vector_space.h"
namespace math_rz
{
	template<typename A,typename B>
	class general_function;

	template<typename A, typename B,typename I=B>
	class integrator;
}
namespace math_rz::analysis::structure::function
{
	template<typename A,math_rz::vector_space_constraint::vector_space  B>
	class norm_topology :public metric_topology<A,B>
	{
	protected:
		using K = typename B::base_field;
		using function_type = math_rz::function<A,B>;
	public:
		virtual real_field norm(const function_type& a) const = 0;
		real_field metric(const function_type& p, const function_type& q) const
		{
			return  norm(p - q);
		}
	};

/*	template<typename A,typename B, int m>
	class Linf_function_norm :public norm_topology<A,B,  m>
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
	template<typename A,typename B, int m>
	class L1_function_norm :public norm_topology<A,B, m>
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
	template<typename A,typename B, int m>
	class Lp_function_norm :public norm_topology<A,B, m>
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

	template<typename A,typename B>
	class Lp_function_norm :public norm_topology<A,B>
	{
		using function_type = math_rz::function<A,B>;
		using K = typename B::K;
		real_field p;
		std::shared_ptr<integrator<A, real_field>> I_ptr;
	public:
		Lp_function_norm(real_field _p, integrator<A, real_field>* _I_ptr) :p(_p),I_ptr(_I_ptr) {}
		Lp_function_norm(real_field _p, std::shared_ptr<integrator<A, real_field>> _I_ptr) 
			:p(_p), I_ptr(_I_ptr) {}


		real_field norm(const function_type& f) const
		{
			return std::pow(I_ptr->integrate
			(
				general_function<A, real_field>([&](const K& u)->K
					{
						return std::pow(f(u).abs(), p);
					})
			), 1. / p);;
		}
	};
}