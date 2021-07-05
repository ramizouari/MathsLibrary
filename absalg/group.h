#pragma once
//#include "maths_concepts.h"
namespace math_rz
{
	class group
	{
	public:
		inline static constexpr enum class Identity0 {} _0 = {};
		group(Identity0)
		{
			//	*this = I0();
		}
		virtual ~group() {};
		/*virtual bool operator==(const group& e) const
		{
			return true;
		}
		bool operator!=(const group& e) const
		{
			return !((*this) == e);
		}*/
		//virtual const group& I0() const =0;

		group& operator+=(Identity0)
		{
			return *this;
		}
		group& operator-=(Identity0)
		{
			return *this;
		}
		virtual bool is_zero() const = 0;
		//virtual group&& operator-() = 0;
		//virtual group& operator==(Identity0) = 0;
	protected:
		constexpr group() {};
	};

}