#pragma once
#include "integrator.h"
namespace math_rz
{
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