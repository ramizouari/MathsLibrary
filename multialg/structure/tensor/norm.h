#pragma once
#include "metric.h"
namespace math_rz::linalg::structure::vector
{
	template<typename K,int n0,int ...n>
	class norm_topology :public metric_topology<K,n>
	{
	protected:
		using tensor_type = metric_topology<K, n0,n...>::vector_type;
	public:
		virtual real_field norm(const tensor_type& a) const = 0;
		real_field metric(const tensor_type& p, const tensor_type& q) const
		{
			return  norm(p - q);
		}
	};

	template<typename K, int n0,...n>
	class Lp_vector_norm :public norm_topology<K,n>
	{
		using tensor_type = norm_topology<K, n0,n...>::vector_type;
		real_field p;
	public:
		Lp_vector_norm(real_field _p) :p(_p) {}

		real_field norm(const tensor_type& M) const
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