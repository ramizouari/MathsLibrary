#pragma once
#include "function.h"
#include "integer.h"

namespace math_rz::analysis {
	template<typename F>
	class sequence :public summable_function<integer, F>
	{
	public:
		sequence(int s0=0) :n0(s0) {}
		F evaluate(const integer& n) const
		{
			return this->operator()(n);
		}

		F partial_sum(const integer& a,const integer& b) const
		{
			F r = 0;
			for (int i = a; i <= b; i++)
				r += evaluate(i);
		}

		F partial_sum(const integer& n) const
		{
			return partial_sum(n0, n);
		}
	private:
		int n0;
	};
}