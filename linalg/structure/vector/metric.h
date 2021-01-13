#pragma once
#include "real_field.h"

namespace math_rz::linalg
{
	template<typename K,int n>
	class finite_dimensional_vector_space;
	template<typename K, int n>
	class square_matrix;
	template<typename K, int n,int m>
	class matrix;
}
namespace math_rz::linalg::structure::vector
{
	template<typename K,int n>
	class metric_topology
	{
	protected:
		using vector_type =math_rz::linalg::finite_dimensional_vector_space<K, n>;
	public:
		virtual real_field metric(const vector_type& p, const vector_type& q) const = 0;
		real_field distance(const vector_type& p, const vector_type& q) const
		{
			return metric(p, q);
		}
	};

	template<typename K,int n>
	class discrete_metric :public metric_topology<K,n>
	{
		using vector_type = metric_topology<K, n>::vector_type;
	public:
		real_field metric(const vector_type& p, const vector_type& q) const
		{
			return p == q ? 1 : 0;
		}
	};

	template<typename K,int n>
	class hamming_metric :public metric_topology<K,n>
	{
		using vector_type = metric_topology<K, n>::vector_type;
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