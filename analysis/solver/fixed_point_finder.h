#pragma once

namespace math_rz::analysis
{
	template<typename E>
	class fixed_point_finder
	{
	public:
		virtual E fixed_point(const function<E, E>& f) const = 0;
	};
}