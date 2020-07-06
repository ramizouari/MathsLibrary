#pragma once
#include "free_algebra.h"
#include "integral_ring.h"
template<typename R>
class polynomial:virtual public free_algebra<R>,virtual public integral_ring
{
public:
	polynomial() {};
	polynomial(const free_algebra<R> &p):free_algebra<R>(p){}
	polynomial(const R &p):free_algebra(p){}
	static polynomial _0()
	{
		return free_algebra<R>::_0();
	}
	static polynomial _1()
	{
		return free_algebra<R>::_1();
	}
	polynomial&& div(const integral_ring& a) { return polynomial(); }
	polynomial&& mod(const integral_ring& b) { return polynomial(); }
	const free_algebra<R>& I0() const
	{
		return free_algebra<R>::I0();
	}
	const free_algebra<R>& I1() const
	{
		return free_algebra<R>::I1();
	}
};

