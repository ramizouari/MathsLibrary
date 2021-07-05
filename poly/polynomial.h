#pragma once
#include <cassert>
#include "free_algebra.h"
#include "absalg/integral_ring.h"
#include "absalg/rational_extension.h"
namespace math_rz::poly {
	namespace structure
	{

		template<typename K>
		class norm_topology;
		template<typename K>
		class inner_product_topology;
	}
	template<typename K>
	class polynomial :virtual public free_algebra<K>, virtual public integral_ring
	{
	public:
		using free_algebra<K>::a;
		using free_algebra<K>::structure_ptr;
		polynomial() {};
		polynomial(const free_algebra<K>& p) :free_algebra<K>(p) {}
		polynomial(free_algebra<K>&& p) :free_algebra<K>(std::move(p)) {}
		polynomial(const std::vector<K>& p) :free_algebra<K>(p) {}
		polynomial(std::vector<K>&& p) :free_algebra<K>(std::move(p)) {}
		polynomial(const K& p) :free_algebra<K>(p) {}
		polynomial(int c) :free_algebra<K>(c) {}
		using base_field = K;
		bool operator!=(const polynomial& p) const
		{
			if (p.degree() != this->degree())
				return true;
			for (int i = 0; i <= p.degree(); i++)
				if (p.a.at(i) != this->a.at(i))
					return true;
			return false;
		}
		bool operator==(const polynomial& p) const
		{
			return !(*this != p);
		}
		static polynomial _0()
		{
			return free_algebra<K>::_0();
		}
		static polynomial _1()
		{
			return free_algebra<K>::_1();
		}

		static std::pair<polynomial, polynomial> euclidean_division(const polynomial& p, const polynomial& q)
		{
			if (p.degree() < q.degree())
				return std::make_pair(0, p);
			polynomial r(p);
			int m(r.degree()), n(q.degree());
			polynomial s;
			s.a.resize(m - n + 1);
			for (; m >= n; m--)
			{
				K k(r.a.at(m) / q.a.at(n));
				s.a.at(m - n) = k;
				if (k.is_zero())
				{
					r.a.pop_back();
					continue;
				}
				for (int i = 1; i <= n; i++)
					r.a.at(m - i) -= k * q.a.at(n - i);
				r.a.pop_back();
			}
			r.reduce();

			return std::make_pair(s, r);
		}
		static polynomial gcd(const polynomial& p, const polynomial& q)
		{
			if (p.degree() < q.degree())
				return gcd(q, p);
			std::pair<polynomial, polynomial> R(euclidean_division(p, q));
			if (R.second.a.empty())
				return q.normalize();
			else return gcd(q, R.second);
		}
		polynomial div(const polynomial& q) const
		{
			return euclidean_division(*this, q).first;
		}
		polynomial mod(const polynomial& q) const
		{
			return euclidean_division(*this, q).second;
		}

		polynomial operator-() const
		{
			return free_algebra<K>::operator-();
		}


		polynomial& operator+=(const polynomial& p)
		{
			this->free_algebra<K>::operator+=(p);
			return *this;
		}
		polynomial& operator-=(const polynomial& p)
		{
			this->free_algebra<K>::operator-=(p);
			return *this;
		}
		polynomial& operator*=(const polynomial& p)
		{
			this->free_algebra<K>::operator*=(p);
			return *this;
		}
		polynomial& operator+=(const K& p)
		{
			this->free_algebra<K>::operator+=(p);
			return *this;
		}
		polynomial& operator-=(const K& p)
		{
			this->free_algebra<K>::operator-=(p);
			return *this;
		}
		polynomial& operator*=(const K& p)
		{
			this->free_algebra<K>::operator*=(p);
			return *this;
		}
		polynomial& operator/=(const K& p)
		{
			for (auto& s : a)
				s /= p;
			this->reduce();
			return *this;
		}

		polynomial& operator/=(int p)
		{
			for (auto& s : a)
				s /= p;
			this->reduce();
			return *this;
		}
		polynomial normalize() const
		{
			polynomial p(*this);
			auto dominant_coeff = p.a.at(p.degree());
			p /= dominant_coeff;
			return p;
		}
		polynomial derivative() const
		{
			polynomial p;
			int n = this->degree();
			p.a.resize(n);
			for (int i = 1; i <= n; i++)
				p.a[i - 1] = this->a[i]*K(i);
			return p;
		}

		polynomial conj() const
		{
			polynomial p = *this;
			for (auto& v : p.a)
				v = v.conj();
			return p;
		}

		real_field norm() const
		{
			return dynamic_cast<poly::structure::norm_topology<K>*>
				(structure_ptr.get())->norm(*this);
		}

		K inner_product(const polynomial& q) const
		{
			return dynamic_cast<poly::structure::inner_product_topology<K>*>
				(structure_ptr.get())->inner_product(*this, q);
		}

		K dot_product(const polynomial& q) const
		{
			return dynamic_cast<poly::structure::inner_product_topology<K>*>
				(structure_ptr.get())->dot_product(*this, q);
		}
		
	};

	template<typename K>
	polynomial<K> operator+(const polynomial<K>& a, const polynomial<K>& b)
	{
		polynomial<K> p(a);
		return p += b;
	}


	template<typename K>
	polynomial<K> operator-(const polynomial<K>& a, const polynomial<K>& b)
	{
		polynomial<K> p(a);
		return p -= b;
	}

	template<typename K>
	polynomial<K> operator*(const polynomial<K>& a, const polynomial<K>& b)
	{
		polynomial<K> p(a);
		return p *= b;
	}

	template<typename K>
	polynomial<K> operator+(const polynomial<K>& a, const K& b)
	{
		polynomial<K> p(a);
		return p += b;
	}

	template<typename K>
	polynomial<K> operator/(const polynomial<K>& a, const K& b)
	{
		polynomial<K> p(a);
		return p /= b;
	}


	template<typename K>
	polynomial<K> operator-(const polynomial<K>& a, const K& b)
	{
		polynomial<K> p(a);
		return p -= b;
	}

	template<typename K>
	polynomial<K> operator*(const polynomial<K>& a, const K& b)
	{
		polynomial<K> p(a);
		return p *= b;
	}

	template<typename K>
	polynomial<K> operator+(const K& b, const polynomial<K>& a)
	{
		polynomial<K> p(a);
		return p += b;
	}


	template<typename K>
	polynomial<K> operator-(const K& b, const polynomial<K>& a)
	{
		polynomial<K> p(a);
		return p -= b;
	}

	template<typename K>
	polynomial<K> operator*(const K& b, const polynomial<K>& a)
	{
		polynomial<K> p(a);
		return p *= b;
	}

	template<typename K>
	polynomial<K> operator/(const polynomial<K>& a, int b)
	{
		polynomial<K> p(a);
		return p *= b;
	}

	template <typename K>
	using rational_function = rational_extension<polynomial<K>>;

	template<typename K>
	std::pair<polynomial<K>, polynomial<K>> bezout(const polynomial<K>& a, const polynomial<K>& b)
	{
		std::pair<polynomial<K>, polynomial<K>> P;
		if (a < b)
		{
			P = bezout(b, a);
			return { P.second,P.first };
		}
		polynomial<K> r0 = a, r1 = b, t0 = 0, t1 = 1, s0 = 1, s1 = 0, w1, w2, w3, q;
		while (!r1.is_zero())
		{
			w1 = r0;
			w2 = t0;
			w3 = s0;
			r0 = r1;
			s0 = s1;
			t0 = t1;
			q = w1.div(r1);
			r1 = w1 - q * r1;
			t1 = w2 - q * t1;
			s1 = w3 - q * s1;
		}
		int n = r0.degree();
		return { s0/r0.coeff(n),t0/r0.coeff(n) };
	}

	template<typename K>
	polynomial<K> gcd(const polynomial<K>& a, const polynomial<K>& b)
	{
		std::pair<polynomial<K>, polynomial<K>> P;
		if (a < b)
		{
			P = bezout(b, a);
			return { P.second,P.first };
		}
		polynomial<K> r0 = a, r1 = b, w1, q;
		while (!r1.is_zero())
		{
			w1 = r0;
			r0 = r1;
			q = w1.div(r1);
			r1 = w1 - q * r1;
		}
		return r0.normalize();
	}
}