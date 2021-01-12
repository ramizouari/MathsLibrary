#pragma once
#include "complex.h"
#include "polynomial.h"

namespace math_rz
{
	template<typename K>
	K root(const polynomial<K>& p,const K& x0,const real_field &err=1e-5)
	{
		auto q = p.derivative();
		K x = x0;
		while (p(x).norm() > err)
			x = x - p(x) / q(x);
		return x;
	}

	template<typename K>
	std::pair<K,integer> root_multiplicity(const polynomial<K>& p, const K& x0, const real_field& err = 1e-5)
	{
		auto w = root(p, x0, err);
		polynomial<K> q = p(polynomial<K>({ -w,1 }));
		int m = 0;
		int n = q.degree();
		while (q.coeff(m).norm()<err)
			m++;
		return { w,m };
	}

	template<typename K>
	std::vector<K> factorise(const polynomial<K>& p, const K& x0, const real_field& eps = 1e-5)
	{
		std::vector<K> f;
		polynomial<K> h = p;
		while (h.degree())
		{
			auto w = root(h, x0, eps);
			h=h.div( polynomial<K>({ -w,1 }));
			f.push_back(w);
		}
		return f;
	}
}