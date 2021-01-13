#pragma once
#include "integrator.h"
namespace math_rz::analysis
{
	/*
	* This class is the base class of "special" integrals
	* - A special integral is an integral which depends on another integral:
	* For example if we work on the euclidean plane, Let's suppose we have an integrator which to every (integrable) function
	* gives its (double) integral over the region R = [0,1]x[0,2pi]
	* To calculate the integral of f over the unit sphere S, We can just calculate the integral of
	* r*f(r*cos(w),r*sin(w)) over the region R
	* We can call the integrator which gives the integral over S of a function f a special integrator depending
	* on the more "elementary" rectangular integrator (or trapezoidal...)
	*/

	/*
	* This class represents integrators depending each one other more simple integrators , see the file
	* for more documentation
	*/
	template<typename E1, typename F1,typename E2,typename F2,typename I=F1>
	class special_integrator : public integrator<E1,F1,I>
	{
	protected:
		std::shared_ptr<integrator<E2, F2,I>> I_ptr;
	public:
		special_integrator(integrator<E2, F2,I>* _I) :I_ptr(_I)
		{

		}

		special_integrator(std::shared_ptr < integrator < E2, F2,I >> _I) :I_ptr(_I)
		{

		}

		void set_integrator(integrator<E2, F2,I>* _I_ptr)
		{
			I_ptr.reset(_I_ptr);
		}

		void set_integrator(std::shared_ptr<integrator<E2, F2,I>> _I_ptr)
		{
			I_ptr.reset(_I_ptr);
		}
	};
}