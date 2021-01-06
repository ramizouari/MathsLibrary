#pragma once
#include "complex.h"
#include "polynomial.h"

namespace math_rz
{
	template<typename F>
	F root(const polynomial<F>& p,const F& x0,const real_field &err=1e-5)
	{
		auto q = p.derivative();
		F x = x0;
		while (p(x).norm() > err)
			x = x - p(x) / q(x);
		return x;
	}

	template<typename F>
	std::pair<F,integer> root_multiplicity(const polynomial<F>& p, const F& x0, const real_field& err = 1e-5)
	{
		auto w = root(p, x0, err);
		polynomial<F> q = p(polynomial<F>({ -w,1 }));
		int m = 0;
		int n = q.degree();
		while (q.coeff(m).norm()<err)
			m++;
		return { w,m };
	}

	template<typename F>
	std::vector<F> factorise(const polynomial<F>& p, const F& x0, const real_field& eps = 1e-5)
	{
		std::vector<F> f;
		polynomial<F> h = p;
		while (h.degree())
		{
			auto w = root(h, x0, eps);
			h=h.div( polynomial<F>({ -w,1 }));
			f.push_back(w);
		}
		return f;
	}
}