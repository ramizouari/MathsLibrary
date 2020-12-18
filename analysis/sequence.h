#pragma once
#include "function.h"
#include "integer.h"

template<typename F>
class sequence :public summable_function<integer,F>
{
public:
	sequence(int s0):n0(s0){}
	F evaluate(int n) const
	{
		return this->operator()(n);
	}
	F partial_sum(int a,int b) const
	{
		F r = 0;
		for (int i = a; i <= b; i++)
			r += evaluate(i);
	}

	F partial_sum(int n) const
	{
		return partial_sum(n0, n);
	}
private:
	int n0;
};