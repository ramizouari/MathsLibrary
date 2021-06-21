#pragma once
#include "real_field.h"

namespace math_rz::multialg
{
	template<typename K,int n0,int ...n>
	class tensor;
}
namespace math_rz::multialg::structure::vector
{
	template<typename K,int n0,int ...n>
	class metric_topology
	{
	protected:
		using tensor_type =tensor<K, n0,n...>;
	public:
		virtual real_field metric(const tensor_type& p, const tensor_type& q) const = 0;
		real_field distance(const tensor_type& p, const tensor_type& q) const
		{
			return metric(p, q);
		}
	};

	template<typename K,int n0,int ...n>
	class discrete_metric :public metric_topology<K,n0,n...>
	{
		using tensor_type = metric_topology<K, n0,n...>::tensor_type;
	public:
		real_field metric(const tensor_type& p, const tensor_type& q) const
		{
			return p == q ? 1 : 0;
		}
	};

	template<typename K,int n0,int ...n>
	class hamming_metric :public metric_topology<K,n0,n...>
	{
		using tensor_type = metric_topology<K, n0,n...> ::vector_type;
	public:
		real_field metric(const tensor_type& p, const tensor_type& q) const
		{
			int d = 0;
			for (int i = 0; i <n; i++)
				if (p[i]!=q[i]) 
					d++;
			return d;
		}
	};

}