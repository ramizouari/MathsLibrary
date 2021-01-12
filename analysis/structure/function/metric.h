#pragma once
#include "real_field.h"

namespace math_rz
{
	template<typename A,typename B>
	class function;
}
namespace math_rz::analysis::structure::function
{
	template<typename A,typename B>
	class metric_topology
	{
	protected:
		using function_type = math_rz::function<A, B>;
	public:
		virtual real_field metric(const function_type& p, const function_type& q) const = 0;
		real_field distance(const function_type& p, const function_type& q) const
		{
			return metric(p, q);
		}
	};

	template<typename A,typename B>
	class discrete_metric :public metric_topology<A,B>
	{
		using function_type = metric_topology<A,B>::function_type;
	public:
		real_field metric(const function_type& p, const function_type& q) const
		{
			return p == q ? 1 : 0;
		}
	};
}