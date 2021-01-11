#pragma once
#include "real_field.h"
#include "linalg/matrix.h"
#include "metric.h"
#include "analysis/integrator/integrator.h"
#include "linalg/structure/vector/norm.h"
namespace math_rz::linalg::structure::matrix
{
	template<typename F,int n,int m>
	class norm_topology :public metric_topology<F,n,m>
	{
	protected:
		using matrix_type = metric_topology<F, n, m>::matrix_type;
	public:
		virtual real_field norm(const matrix_type& a) const = 0;
		real_field metric(const matrix_type& p, const matrix_type& q) const
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

	template<typename F, int n, int m>
	class Lpq_vector_norm :public norm_topology<F,n,m>
	{
		using matrix_type = norm_topology<F, n, m>::matrix_type;
		real_field p,q;
	public:
		Lpq_vector_norm(real_field _p, real_field _q) :p(_p),q(_q) {}

		real_field norm(const matrix_type& M) const
		{
			real_field R = 0,r;
			constexpr auto inf = std::numeric_limits<long double>::infinity();
			for (int i = 0; i < n; i++) 
			{
				r = 0;
				for (int j = 0; j < m; j++)
					if (q < inf)
						r += std::pow(M[i][j].abs(), q);
					else r = std::max(r, M[i][j]);
				if (q < inf)
				{
					if (p < inf)
						R += std::pow(r, p / q);
					else R = std::max<long double>(R, std::pow(r, 1. / q));
				}
				else
				{
					if (p < inf)
						R += std::pow(r,p);
					else R = std::max(R, r);
				}
			}
			if (p < inf)
				R = std::pow(R, 1. / p);
			return R;
		}
	};

	template<typename F, int n, int m>
	class L22_operator_norm :public norm_topology<F, n, m>
	{
		using matrix_type = norm_topology<F, n, m>::matrix_type;
	public:
		L22_operator_norm()  {}

		real_field norm(const matrix_type& M) const
		{
			return largest_sing(M);
		}
	};

	template<typename F, int n, int m>
	class L1q_operator_norm :public norm_topology<F, n, m>
	{
		using matrix_type = norm_topology<F, n, m>::matrix_type;
		real_field q;
	public:
		L1q_operator_norm(real_field _q):q(_q) {}

		real_field norm(const matrix_type& M) const
		{
			math_rz::linalg::structure::vector::Lp_vector_norm<F,n> N(q);
			auto M_t = M.transpose();
			F result;
			for (int i = 0; i < m; i++)
			{
				coordinate_space<F, n>&& u(std::move(M_t[i]));
				result = std::max(result, N.norm(u));
			}
			return result;
		}
	};

	template<typename F,int n,int m>
	class Lpinf_operator_norm :public norm_topology<F, n, m>
	{
		using matrix_type = norm_topology<F, n, m>::matrix_type;
		real_field p;
	public:
		Lpinf_operator_norm(real_field _p) :p(_p) {}

		real_field norm(const matrix_type& M) const
		{
			constexpr auto inf = std::numeric_limits<long double>::infinity();
			static L1q_operator_norm<F, n, m> N1_inf(inf);
			static math_rz::linalg::structure::vector::Lp_vector_norm<F, n> N(inf);
			if (p == inf)
				return N1_inf.norm(M);
			F result;
			for (int i = 0; i < m; i++)
			{
				coordinate_space<F, n> u(M[i]);
				result = std::max(result, N.norm(u));
			}
			return result;
		}
	};
}