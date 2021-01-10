#pragma once
#include "absalg/ring.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <execution>
#include <compare>
#include "real_field.h"

namespace math_rz {
	namespace poly::multiplicator 
	{
		template<typename F>
		class multiplicator;

		template<typename F>
		class karatsuba_multiplicator;
	}

	namespace poly::structure
	{
		template<typename F>
		class metric_topology;
	}

	template<typename R>
	class free_algebra : virtual public ring
	{
	public:
		free_algebra() {}
		free_algebra(R m) :a(1, m) { reduce(); }
		free_algebra(std::vector<R>&& c) :a(std::move(c)) { reduce(); }
		free_algebra(const std::vector<R>& c) :a(c) { reduce(); }
		free_algebra(const free_algebra<R>& p) :a(p.a) {}
		free_algebra(int c) :a(1, R(c)) {}


		bool operator<(const free_algebra& q) const
		{
			return degree() < q.degree();
		}

		int degree() const
		{
			return a.size() - 1;
		}

		free_algebra operator-() const
		{
			free_algebra p(*this);
			for (auto& s : p.a)
				s = -s;
			return p;
		}

		template<typename H = R>
		H operator()(const H& u) const
		{
			H r = 0, w = 1;
			for (const auto& s : a)
			{
				r += w * s;
				w *= u;
			}
			return r;
		}

		const free_algebra& operator+() const
		{
			return *this;
		}

		free_algebra& operator+=(const free_algebra& p)
		{
			for (int i = 0; i <= degree(); i++)
				if (i > p.degree())
					break;
				else a[i] += p.a[i];
			for (int i = degree() + 1; i <= p.degree(); i++)
				a.push_back(p.a[i]);
			reduce();
			return *this;
		}
		free_algebra& operator+=(const R& p)
		{
			return *this += free_algebra(p);
		}
		free_algebra& operator-=(const R& p)
		{
			return *this += free_algebra(p);
		}

		free_algebra& operator*=(const R& p)
		{
			std::for_each(a.begin(), a.end(), [&p](auto& v) {v *= p; });
			reduce();
			return *this;
		}
		free_algebra& operator-=(const free_algebra& p)
		{
			for (int i = 0; i <= degree(); i++)
				if (i > p.degree())
					break;
				else a[i] -= p.a[i];
			for (int i = degree() + 1; i <= p.degree(); i++)
				a.push_back(p.a[i]);
			reduce();
			return *this;
		}
		free_algebra& operator*=(const free_algebra& p)
		{
			return *this = multiply_ptr->multiply(*this, p);
		}
		const R& coeff(int n) const
		{
			return a.at(n);
		}

		R& coeff(int n)
		{
			return a.at(n);
		}

		virtual ring& operator+=(int n)
		{
			return *this += R(n);
		}
		virtual ring& operator-=(int n)
		{
			return *this -= R(n);
		}
		virtual ring& operator*=(int n)
		{
			return *this *= R(n);
		}
		bool is_zero() const
		{
			//assert(!(a.size() == 1 && a.back().is_zero()));
			return a.empty();
		}
		bool is_one() const
		{
			return !a.empty() && a.at(0).is_one();
		}

		std::vector<R>& get_vect()
		{
			return a;
		}

		const std::vector<R>& get_vect() const
		{
			return a;
		}
		void reduce() {
			while (!a.empty() && (a.back().is_zero()))
				a.pop_back();
		}
		
		static const math_rz::poly::multiplicator::multiplicator<R>& get_multiplicator()
		{
			return *multiply_ptr;
		}

		static void set_multiplicator(math_rz::poly::multiplicator::multiplicator<R>* m_ptr)
		{
			multiply_ptr.reset(m_ptr);
		}


		static void set_structure(poly::structure::metric_topology<R>* N)
		{
			structure_ptr.reset(N);
		}

		static const poly::structure::metric_topology<R>& get_structure()
		{
			return *structure_ptr;
		}

		real_field metric(const free_algebra& p)
		{
			return structure_ptr->metric(*this, p);
		}

		real_field distance(const free_algebra& p)
		{
			return metric(p);
		}

		protected:
		std::vector<R> a;

		inline static std::unique_ptr<math_rz::poly::multiplicator::multiplicator<R>> 
			multiply_ptr =std::make_unique< math_rz::poly::multiplicator::karatsuba_multiplicator<R>>
			(math_rz::poly::multiplicator::karatsuba_multiplicator<R>());

		inline static std::unique_ptr<poly::structure::metric_topology<R>> structure_ptr;
	};

	template<typename R>
	free_algebra<R> operator+(const free_algebra<R>& a, const free_algebra<R>& b)
	{
		free_algebra<R> p(a);
		return p += b;
	}


	template<typename R>
	free_algebra<R> operator-(const free_algebra<R>& a, const free_algebra<R>& b)
	{
		free_algebra<R> p(a);
		return p -= b;
	}

	template<typename R>
	free_algebra<R> operator*(const free_algebra<R>& a, const free_algebra<R>& b)
	{
		free_algebra<R> p(a);
		return p *= b;
	}

	template<typename R>
	free_algebra<R> operator+(const free_algebra<R>& a, const R& b)
	{
		free_algebra<R> p(a);
		return p += b;
	}


	template<typename R>
	free_algebra<R> operator-(const free_algebra<R>& a, const R& b)
	{
		free_algebra<R> p(a);
		return p -= b;
	}

	template<typename R>
	free_algebra<R> operator*(const free_algebra<R>& a, const R& b)
	{
		free_algebra<R> p(a);
		return p *= b;
	}

	template<typename R>
	free_algebra<R> operator+(const R& b, const free_algebra<R>& a)
	{
		free_algebra<R> p(a);
		return p += b;
	}


	template<typename R>
	free_algebra<R> operator-(const R& b, const free_algebra<R>& a)
	{
		free_algebra<R> p(a);
		return p -= b;
	}

	template<typename R>
	free_algebra<R> operator*(const R& b, const free_algebra<R>& a)
	{
		free_algebra<R> p(a);
		return p *= b;
	}

	template<typename R>
	std::ostream& operator<<(std::ostream& H, const free_algebra<R>& p)
	{
		H << "(";
		for (int i = 0; i <= p.degree(); i++)
			if (i < p.degree())
				H << p.coeff(i) << ", ";
			else H << p.coeff(i) << ")";
		return H;
	}
}