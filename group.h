#pragma once
//#include "maths_concepts.h"
class group
{
	public:
	static enum class Identity0 {} _0;
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
	//virtual group&& operator-() = 0;
	//virtual group& operator==(Identity0) = 0;
protected:
	constexpr group() {};
};

