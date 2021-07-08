#pragma once
//#include "maths_concepts.h"
namespace math_rz
{
	class group
	{
	public:
		virtual ~group() {};

		virtual bool is_zero() const = 0;
		//virtual group&& operator-() = 0;
		//virtual group& operator==(Identity0) = 0;
	protected:
		constexpr group() {};
	};

}