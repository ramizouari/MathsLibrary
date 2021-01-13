#pragma once

namespace math_rz::analysis
{
	template<typename A>
	A exp(const A& w,const integer& n)
	{
		A r,u(1);
		for (int i = 0; i <= n; i++, u *= w / A(i))
			r += u;
		return r;
	}

	template<typename A>
	A exp(const A& w, long long n)
	{
		return exp(w, integer(n));
	}

	template<typename A>
	A exp(const A& w, int n)
	{
		return exp(w, integer(n));
	}

	template<typename A>
	A exp_m1(const A& w, const integer& n)
	{
		A r, u(1);
		int i = 0;
		for (int i = 0; i <= n; i++, u *= w / A(i))if (i == 0) continue;
			else r += u;
		return r;
	}

	template<typename A>
	A exp_m1(const A& w,long long n)
	{
		return exp_m1(w, integer(n));
	}

	template<typename A>
	A exp_m1(const A& w, int n)
	{
		return exp_m1(w, integer(n));
	}

	template<typename A>
	A exp(const A& w,const real_field &err)
	{
		A r, u(1);
		int i = 0;
		for (i = 0; u.norm()-r.norm()*err>0; i++, u *= w / A(i))
			r += u;
		return r;
	}

	template<typename A>
	A exp(const A& w, long double err = 1e-5)
	{
		return exp(w, real_field(err));
	}

	template<typename A>
	A exp(const A& w,  double err )
	{
		return exp(w, real_field(err));
	}


	template<typename A>
	A exp_m1(const A& w, const real_field& err)
	{
		A r, u(1);
		int i = 0;
		for (i = 0; u.norm() - r.norm() * err > 0; i++, u *= w / A(i))if (i == 0) continue;
			else r += u;
		return r;
	}

	template<typename A>
	A exp_m1(const A& w, long double err = 1e-5)
	{
		return exp_m1(w, real_field(err));
	}

	template<typename A>
	A exp_m1(const A& w, double err = 1e-5)
	{
		return exp_m1(w, real_field(err));
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
	auto pow(const A& u, const B& v,int n)
	{
		return exp(log(u, n) * v, n);
	}

}