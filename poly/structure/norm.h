#pragma once
#include "real_field.h"
#include "poly/polynomial.h"
#include "metric.h"
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
	class Linf_norm:public norm_topology<F>
	{
		real_field a, b;
		integer cuts;
	public:
		Linf_norm(real_field _a,real_field _b,integer _c):a(_a),b(_b),cuts(_c){}
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
		real_field a, b;
		integer cuts;
	public:
		L1_function_norm(real_field _a, real_field _b, integer _c) :a(_a), b(_b), cuts(_c) {}
		real_field norm(const polynomial<F>& f) const
		{
			real_field eps = (b - a) / cuts;
			real_field u = a, R = 0;
			for (int i = 0; i <= cuts; i++, u += eps)
				R += f(u).abs();
			return R;
		}
	};

	template<typename F>
	class Lp_function_norm :public norm_topology<F>
	{
		real_field a, b;
		integer cuts,p;
	public:
		Lp_function_norm(integer _p,real_field _a, real_field _b, integer _c)
			:a(_a), b(_b), cuts(_c),p(_p) {}
		
		real_field norm(const polynomial<F>& f) const
		{
			real_field eps = (b - a) / cuts;
			real_field u = a, R = 0;
			for (int i = 0; i <= cuts; i++, u += eps)
				R += std::pow(f(u).abs(),p);
			return std::pow(R,1./p);
		}
	};

	template<typename F>
	class Lp_vector_norm :public norm_topology<F>
	{
		integer p;
	public:
		Lp_vector_norm(integer _p):p(_p) {}

		real_field norm(const polynomial<F>& f) const
		{
			real_field R = 0;
			for (const auto& v : f.get_vect())
				R += std::pow(v.abs(), p);
			return std::pow(R, 1. / p);
		}
	};
}