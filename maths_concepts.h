#pragma once
#include <concepts>
namespace maths_concepts
{
	template<typename G>
	concept monoid =  std::equality_comparable<G>&& requires(G a,G b)
	{
		a * b->G;
	};
}