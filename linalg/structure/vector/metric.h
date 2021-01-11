#pragma once
#include "real_field.h"

namespace math_rz
{
	template<typename F,int n>
	class finite_dimensional_vector_space;
	template<typename F, int n>
	class square_matrix;
	template<typename F, int n,int m>
	class matrix;
}
namespace math_rz::linalg::structure::vector
{
	template<typename F,int n>
	class metric_topology
	{
	protected:
		using vector_type =math_rz::finite_dimensional_vector_space<F, n>;
	public:
		virtual real_field metric(const vector_type& p, const vector_type& q) const = 0;
		real_field distance(const vector_type& p, const vector_type& q) const
		{
			return metric(p, q);
		}
	};

	template<typename F,int n>
	class discrete_metric :public metric_topology<F,n>
	{
		using vector_type = metric_topology<F, n>::vector_type;
	public:
		real_field metric(const vector_type& p, const vector_type& q) const
		{
			return p == q ? 1 : 0;
		}
	};

	template<typename F,int n>
	class hamming_metric :public metric_topology<F,n>
	{
		using vector_type = metric_topology<F, n>::vector_type;
	public:
		real_field metric(const vector_type& p, const vector_type& q) const
		{
			int d = 0;
			for (int i = 0; i <n; i++)
				if (p[i]!=q[i]) 
					d++;
			return d;
		}
	};

}