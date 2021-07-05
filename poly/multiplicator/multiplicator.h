#pragma once
#include "poly/polynomial.h"
#include "poly/interpolation.h"
#include <future>
#include <atomic>
namespace math_rz::poly::multiplicator
{
	template<typename R>
	class multiplicator
	{
	public:
		virtual free_algebra<R> multiply(const free_algebra<R>& p,const  free_algebra<R>& q) const
		{
			if (p.degree() < 0 || q.degree() < 0)
				return 0;
			int m = p.degree() + q.degree();
			free_algebra<R> h;
			h.get_vect().resize(m + 1);
			for (int i = 0; i <= p.degree(); i++)
				for (int j = 0; j <= q.degree(); j++)
					h.coeff(i + j) += p.coeff(i) * q.coeff(j);
			h.reduce();
			return h;
		}
	};

	template<typename R>
	class karatsuba_multiplicator:public multiplicator<R>
	{
	public:
		inline static constexpr int limit = 50;
		virtual free_algebra<R> multiply(const free_algebra<R>& p, const free_algebra<R>& q) const
		{
			int n = p.degree(), m = q.degree();
			int v = std::min(n, m), w = std::max(n, m);
			if (w < limit)
				return multiplicator<R>::multiply(p, q);
			int s = w / 2;
			free_algebra<R> p1, p2, q1, q2;
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
			free_algebra<R> d = multiply(p1 + p2, q1 + q2), a = multiply(p1, q1), c = multiply(p2, q2);
			free_algebra<R> b = d - a - c;
			auto k1 = a.degree(), k2 = b.degree(), k3 = c.degree();
			auto K = std::max({ k3 + 2 * s,k2 + s,k1 });
			free_algebra<R> r;
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
	};


	template<typename R>
	class parallel_karatsuba_multiplicator :public karatsuba_multiplicator<R>
	{
		mutable std::atomic<int> threads=0;
		int coeff;
	public:
		parallel_karatsuba_multiplicator(int _coeff=1):coeff(_coeff){}
		inline static constexpr int limit = 50;
		virtual free_algebra<R> multiply(const free_algebra<R>& p, const free_algebra<R>& q) const
		{
			threads++;
			if (threads > coeff * std::thread::hardware_concurrency())
			{
				threads--;
				return karatsuba_multiplicator<R>::multiply(p, q);
			}
			int n = p.degree(), m = q.degree();
			int v = std::min(n, m), w = std::max(n, m);
			if (w < limit)
			{
				threads--;
				return multiplicator<R>::multiply(p, q);
			}
			int s = w / 2;
			free_algebra<R> p1, p2, q1, q2;
			p1.get_vect().resize(std::min(s, n + 1));
			q1.get_vect().resize(std::min(s, m + 1));
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
			free_algebra<R>  d, a, c;
			auto C1 = std::async(&parallel_karatsuba_multiplicator::multiply,this,p1 + p2, q1 + q2), 
				C2 =std::async (&parallel_karatsuba_multiplicator::multiply,this,p1, q1), 
				C3 = std::async(&parallel_karatsuba_multiplicator::multiply,this,p2, q2);
			d = C1.get(); 
			a = C2.get();
			c = C3.get();
			free_algebra<R> b = d - a - c;
			auto k1 = a.degree(), k2 = b.degree(), k3 = c.degree();
			auto K = std::max({ k3 + 2 * s,k2 + s,k1 });
			free_algebra<R> r;
			r.get_vect().resize(K + 1, 0);
			for (int i = 0; i <= k1; i++)
				r.coeff(i) += a.coeff(i);
			for (int i = s; i <= s + k2; i++)
				r.coeff(i) += b.coeff(i - s);
			for (int i = 2 * s; i <= 2 * s + k3; i++)
				r.coeff(i) += c.coeff(i - 2 * s);
			r.reduce();
			threads--;
			return r;
		}
	};

	class fast_multiplicator :public multiplicator<complex>
	{
		using R = complex;
	public:
		
		inline static constexpr int limit = 50;
		virtual free_algebra<R> multiply(const free_algebra<R>& p, const free_algebra<R>& q) const
		{
			
		}
	};
}