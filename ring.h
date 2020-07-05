#pragma once
#include "group.h"
class ring :
    virtual public group
{
public:
	static enum class Identity1 {} _1;
	using Identity0 = group::Identity0;
	virtual ~ring() {};
	ring(Identity0) :group(_0)
	{
	}
	ring(Identity1)
	{
		*this = I1();
	}
	virtual const ring& I1() const = 0;
	virtual const ring& I0() const = 0;
	ring& operator*=(Identity1)
	{
		return *this;
	}
	ring& operator/=(Identity1)
	{
		return *this;
	}
	ring& operator*=(Identity0)
	{
		*this = I0();
		return *this;
	}
protected:
	ring() {};
};

namespace algebra
{
	template<typename R>
	R pow(const R& u, int n)
	{
		if (n == 0)
			return R::_1();
		R v = pow<R>(u, n / 2);
		R w(v * v);
		if ((n % 2) == 0)
			return w;
		else
			return R(w * u);
	}

	
}