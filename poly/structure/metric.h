#pragma once
#include "real_field.h"
#include "poly/polynomial.h"

namespace math_rz::poly::structure
{
	template<typename F>
	class metric_topology
	{
	public:
		virtual real_field metric(const free_algebra<F>& p, const free_algebra<F>& q) const = 0;
		real_field distance(const free_algebra<F>& p, const free_algebra<F>& q) const
		{
			return metric(p, q);
		}
	};

	template<typename F>
	class discrete_metric:public metric_topology<F>
	{
	public:
		real_field metric(const free_algebra<F>& p, const free_algebra<F>& q) const
		{
			return p == q ? 1 : 0;
		}
	};

	template<typename F>
	class hamming_metric :public metric_topology<F>
	{
		real_field metric(const free_algebra<F>& p, const free_algebra<F>& q) const
		{
			int d = 0;
			int n = p.degree(), m = q.degree(),s=std::min(n,m);
			for (int i = 0; i <= s; i++) if (p.coeff(i) != q.coeff(i)) d++;
			return d + (std::max(n, m) - s);
		}
	};
	
}