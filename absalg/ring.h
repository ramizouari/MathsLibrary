#pragma once
#include "group.h"
#include <concepts>
namespace math_rz
{
	class ring :
		virtual public group
	{
	public:
		virtual ~ring() {};
		//virtual const ring& I1() const = 0;
		//virtual const ring& I0() const = 0;
		virtual bool is_one() const = 0;
		//virtual ring& operator+=(int n) = 0;
		//virtual ring& operator-=(int n) = 0;
		//virtual ring& operator*=(int n) = 0;
	protected:
		ring() {};
	};

	namespace ring_constraints
	{

		template<typename G>
		concept group = requires (const G & a, const G & b)
		{
			{a* b}->std::convertible_to<G>;
			{a.inv()}->std::convertible_to<G>;
		};

		template<typename G>
		concept commutative_group = requires (const G & a, const G & b)
		{
			{a + b}->std::convertible_to<G>;
			{a - b}->std::convertible_to<G>;
			{-a}->std::convertible_to<G>;
		};

		template<typename G>
		concept ring = std::is_base_of_v<math_rz::ring, G> && commutative_group<G> && requires (const G & a, const G & b)
		{
			{a* b}->std::convertible_to<G>;
		};
		template<typename G>
		concept field = ring<G> && commutative_group<G> && group<G>;

		template<typename S>
		concept ordered =  requires(const S & a, const S & b)
		{
			{a <=> b};
		};

		template<typename R>
		concept has_abs = requires(R a)
		{
			a.abs();
		};
	}

	template<typename A>
	A unit(const A&u)
	{
		if constexpr (ring_constraints::ring<A>)
			return 1;
		else return A(u, 1);
	}

	template<typename A>
	A unit()
	{
		return 1;
	}

	template<typename A>
	A zero(const A& u)
	{
		if constexpr (ring_constraints::ring<A>)
			return A();
		else return A(u, 0);
	}

	template<typename A>
	A zero()
	{
		return A();
	}

	template<typename A>
	A pow(const A& u, long long n)
	{
		if (n == 0)
			return unit(u);
		else if (n == 1)
			return u;
		A r = pow(u, n / 2);
		r *= r;
		return r * pow(u, n % 2);
	}

	template<typename A>
	A commutator(const A& u,const A& v)
	{
		return u * v - v * u;
	}
}