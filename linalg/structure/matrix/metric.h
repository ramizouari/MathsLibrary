#pragma once
#include "real_field.h"
namespace math_rz::linalg
{
	template<typename K,int n,int m>
	class matrix;
}
namespace math_rz::linalg::structure::matrix
{
	template<typename K,int n,int m>
	class metric_topology
	{
	protected:
	public:
		using matrix_type = math_rz::linalg::matrix<K, n, m>;
		virtual real_field metric(const matrix_type& p, const matrix_type& q) const = 0;
		real_field distance(const matrix_type& p, const matrix_type& q) const
		{
			return metric(p, q);
		}
	};

	template<typename K,int m,int n>
	class discrete_metric :public metric_topology<K,n,m>
	{
	public:
		using matrix_type = math_rz::linalg::matrix<K, n, m>;
		real_field metric(const matrix_type& p, const matrix_type& q) const
		{
			return p == q ? 1 : 0;
		}
	};

	template<typename K,int n,int m>
	class hamming_metric :public metric_topology<K,n,m>
	{
		using matrix_type = metric_topology<K, n, m>::matrix_type;
	public:
		real_field metric(const matrix_type& p, const matrix_type& q) const
		{
			int d = 0;
			for (int i = 0; i <n; i++)
				for(int j=0;j<m;j++)
					if (p[i][j]!=q[i][j]) 
						d++;
			return d;
		}
	};

}