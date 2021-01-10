#pragma once
#include <random>
namespace math_rz
{
	template<typename S>
	class generator
	{
	public:
		virtual S generate() = 0;

	};
}