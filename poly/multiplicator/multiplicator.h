#pragma once
#include "poly/polynomial.h"
namespace math_rz::poly::multiplicator
{
	template<typename F>
	class multiplicator
	{
	public:
		virtual free_algebra<F> multiply(const free_algebra<F>& p,const  free_algebra<F>& q) const
		{
			if (p.degree() < 0 || q.degree() < 0)
				return 0;
			int m = p.degree() + q.degree();
			free_algebra<F> h;
			h.get_vect().resize(m + 1);
			for (int i = 0; i <= p.degree(); i++)
				for (int j = 0; j <= q.degree(); j++)
					h.coeff(i + j) += p.coeff(i) * q.coeff(j);
			h.reduce();
			return h;
		}
	};

	template<typename F>
	class karatsuba_multiplicator:public multiplicator<F>
	{
	public:
		inline static constexpr int limit = 50;
		virtual free_algebra<F> multiply(const free_algebra<F>& p, const free_algebra<F>& q) const
		{
			int n = p.degree(), m = q.degree();
			int v = std::min(n, m), w = std::max(n, m);
			if (w < limit)
				return multiplicator<F>::multiply(p, q);
			int s = w / 2;
			free_algebra<F> p1, p2, q1, q2;
			p1.get_vect().resize(std::min(s,n+1));
			q1.get_vect().resize(std::min(s,m+1));
			p2.get_vect().resize(std::max(n + 1 - s, 0));
			q2.get_vect().resize(std::max(m + 1 - s, 0));
			for (int i = 0; i < std::min(s, n + 1); i++)
				p1.coeff(i) = p.coeff(i);

			for (int i = 0; i < std::min(s, m + 1); i++)
				q1.coeff(i) = q.coeff(i);

			for (int i = s; i <= n; i++)
				p2.coeff(i - s) = p.coeff(i);
			for (int i = s; i <= m; i++)
				q2.coeff(i - s) = q.coeff(i);
			free_algebra<F> d = multiply(p1 + p2, q1 + q2), a = multiply(p1, q1), c = multiply(p2, q2);
			free_algebra<F> b = d - a - c;
			auto k1 = a.degree(), k2 = b.degree(), k3 = c.degree();
			auto K = std::max({ k3 + 2 * s,k2 + s,k1 });
			free_algebra<F> r;
			r.get_vect().resize(K + 1, 0);
			for (int i = 0; i <= k1; i++)
				r.coeff(i) += a.coeff(i);
			for (int i = s; i <= s + k2; i++)
				r.coeff(i) += b.coeff(i - s);
			for (int i = 2 * s; i <= 2 * s + k3; i++)
				r.coeff(i) += c.coeff(i - 2 * s);
			r.reduce();
			return r;
		}

		/*
		* virtual free_algebra<F> multiply(const polynomial<F>& p, const polynomial<F>& q) const
		{
			int n = p.degree(), m = q.degree();
			int w = std::min(n, m);
			if (w < 50)
				return multiplicator<F>::multiply(p, q);
			int s = w / 2;
			free_algebra<F> p1, p2, q1, q2;
			p1.get_vect().resize(s);
			q1.get_vect().resize(s);
			p2.get_vect().resize(n+1 - s);
			q2.get_vect().resize(m + 1 - s);
			for (int i = 0; i < s; i++)
			{
				p1.coeff(i) = p.coeff(i);
				q1.coeff(i) = q.coeff(i);
			}
			for (int i = s; i <= n; i++)
				p2.coeff(i - s) = p.coeff(i);
			for (int i = s; i <= m; i++)
				q2.coeff(i - s) = q.coeff(i);
			free_algebra<F> d = multiply(p1 + p2, q1 + q2),a=multiply(p1,q1),c=multiply(p2,q2);
			free_algebra<F> b = d - a - c;
			auto k1 = a.degree(), k2=b.degree(), k3=c.degree();
			auto K = std::max({ k3 + 2 * s,k2 + s,k1});
			free_algebra<F> r;
			r.get_vect().resize(K+1,0);
			for (int i = 0; i <= k1; i++)
				r.coeff(i) += a.coeff(i);
			for (int i = s; i <= s+k2; i++)
				r.coeff(i) += b.coeff(i - s);
			for (int i = 2*s; i <= 2*s+k3; i++)
				r.coeff(i) += c.coeff(i - 2*s);
			return r;
		}
		*/
	};
}