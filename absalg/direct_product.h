#pragma once
#include "ring.h"
#include <vector>
namespace math_rz::absalg
{
	template<ring_constraints::group G, int n>
	class direct_product
	{
		std::vector<G> u;
	public:
		direct_product(const std::vector<G> &P):u(P){}
		direct_product& operator+=(const direct_product& o)
		{
			for (int i = 0; i < n; i++)
				u[i] += o.u[i];
			return *this;
		}
		direct_product& operator-=(const direct_product& o)
		{
			for (int i = 0; i < n; i++)
				u[i] += o.u[i];
			return *this;
		}
		direct_product& operator*=(const direct_product& o) requires ring_constraints::ring<G>
		{
			for (int i = 0; i < n; i++)
				u[i] += o.u[i];
			return *this;
		}

		direct_product operator+(const direct_product& o)
		{
			direct_product P = *this;
			return P += o;
		}
		direct_product operator-(const direct_product& o)
		{
			direct_product P = *this;
			return P -= o;
		}
		direct_product operator*(const direct_product& o) requires ring_constraints::ring<G>
		{
			direct_product P = *this;
			return P *= o;
		}
	};
}
