#pragma once
#include "metric.h"
namespace math_rz::linalg::structure::vector
{
	template<typename K,int n>
	class norm_topology :public metric_topology<K,n>
	{
	protected:
		using vector_type = metric_topology<K, n>::vector_type;
	public:
		virtual real_field norm(const vector_type& a) const = 0;
		real_field metric(const vector_type& p, const vector_type& q) const
		{
			return  norm(p - q);
		}
	};

/*	template<typename K, int n, int m>
	class Linf_function_norm :public norm_topology<K, n,  m>
	{
		real_field a, b;
		integer cuts;
		
	public:
		Linf_norm(real_field _a, real_field _b, integer _c) :a(_a), b(_b), cuts(_c) {}
		real_field norm(const polynomial<K>& f) const
		{
			real_field eps = (b - a) / cuts;
			real_field u = a, R = 0;
			for (int i = 0; i <= cuts; i++, u += eps)
				R = std::max({ R,f(u).abs() });
			return R;
		}
	};*/
/*
	template<typename K, int n, int m>
	class L1_function_norm :public norm_topology<K, n, m>
	{
		std::shared_ptr<integrator<K, K>> I_ptr;
	public:
		L1_function_norm(std::shared_ptr<integrator<K, K>> _I_ptr) :I_ptr(_I_ptr) {}
		L1_function_norm(integrator<K, K>* _I_ptr) :I_ptr(_I_ptr) {}
		real_field norm(const polynomial<K>& f) const
		{
			return I_ptr->integrate
			(
				general_function<K, K>([&](const K& u)->K
					{
						return f(u).abs();
					})
			);
		}
	};*/
	/*
	template<typename K, int n, int m>
	class Lp_function_norm :public norm_topology<K, n, m>
	{
		integer p;
		std::shared_ptr<integrator<K, K>> I_ptr;
	public:
		Lp_function_norm(integer _p, std::shared_ptr<integrator<K, K>> _I_ptr) :I_ptr(_I_ptr), p(_p) {}
		Lp_function_norm(integer _p, integrator<K, K>* _I_ptr) :I_ptr(_I_ptr), p(_p) {}

		real_field norm(const polynomial<K>& f) const
		{
			return std::pow(I_ptr->integrate
			(
				general_function<K, K>([&](const K& u)->K
					{
						return std::pow(f(u).abs(), p);
					})
			), 1. / p);
		}
	};*/

	template<typename K, int n>
	class Lp_vector_norm :public norm_topology<K,n>
	{
		using vector_type = norm_topology<K, n>::vector_type;
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
				else R = std::max(R, M[i].abs());
			if (p < inf)
				R = std::pow(R, 1. / p);
			return R;
		}
	};
}