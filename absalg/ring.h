#pragma once
#include "group.h"
#include <concepts>
namespace math_rz
{
	class ring :
		virtual public group
	{
	public:
		static enum class Identity1 {} _1;
		using Identity0 = group::Identity0;
		virtual ~ring() {};
		ring(Identity0) :group(_0)
		{
		}
		ring(Identity1)
		{
			//	*this = I1();
		}
		//virtual const ring& I1() const = 0;
		//virtual const ring& I0() const = 0;
		virtual bool is_one() const = 0;
		//virtual ring& operator+=(int n) = 0;
		//virtual ring& operator-=(int n) = 0;
		//virtual ring& operator*=(int n) = 0;
	protected:
		ring() {};
	};


	template<typename A>
	A pow(const A& u, long long n)
	{
		if (n == 0)
			return 1;
		else if (n == 1)
			return u;
		//else if (n < 0)
		//	return A(1) / pow(u, -n);
		A r = pow(u, n / 2);
		r *= r;
		return r * pow(u, n % 2);
	}

	

	template<typename A>
	A commutator(const A& u,const A& v)
	{
		return u * v - v * u;
	}


	namespace ring_constraints
	{

		template<typename G>
		concept group = requires (const G & a, const G & b)
		{
			{a* b}->std::convertible_to<G>;
			{inv(a)}->std::convertible_to<G>;
		};

		template<typename G>
		concept commutative_group = requires (const G & a, const G & b)
		{
			{a + b}->std::convertible_to<G>;
			{a - b}->std::convertible_to<G>;
			{-a}->std::convertible_to<G>;
		};

		template<typename G>
		concept ring = commutative_group<G> && requires (const G & a, const G & b)
		{
			{a * b}->std::convertible_to<G>;
		};
		template<typename G>
		concept field = ring<G> && commutative_group<G> && group<G>;

		template<typename S>
		concept ordered = requires(const S & a, const S & b)
		{
			{a <=> b};
		};
	
		template<typename R>
		concept has_abs = requires(R a)
		{
			a.abs();
		};
	}
}