#pragma once

namespace math_rz
{
	template<typename A>
	A exp(const A& w,int n=8)
	{
		A r,u(1);
		for (int i = 0; i <= n; i++, u *= w / A(i))
			r += u;
		return r;
	}

	template<typename A>
	A geom_sum(const A& w, int n)
	{
		A r, u(1);
		for (int i = 0; i <= n; i++, u *= w)
			r += u;
		return r;
	}

	template<typename A>
	A geom_sum(const A& w)
	{
		return A(1)/(A(1) - w);
	}


	template<typename A>
	A log(const A& w, int n)
	{
		A r, u(w-A(1)),v(w-A(1));
		for (int i = 1; i <= n; i++, u *= v)
			if (i % 2 == 0)
				r -= u / A(i);
			else r += u/A(i);
		return r;
	}

	template<typename A,typename B>
	auto pow(const A& u, const B& v,int n=8)
	{
		return math_rz::exp(math_rz::log(u, n) * v, n);
	}

	template<typename A>
	A pow(const A& u, int n)
	{
		if (n == 0)
			return 1;
		else if (n == 1)
			return u;
		else if (n < 0)
			return A(1) / pow(u, -n);
		A r = pow(u, n / 2);
		r *= r;
		return r * pow(u, n % 2);
	}
}