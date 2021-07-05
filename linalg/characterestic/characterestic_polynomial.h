#pragma once
#include "linalg/matrix.h"
#include "poly/polynomial.h"
#include "absalg/ring.h"
namespace math_rz::linalg::characterestic
{
	template<typename R,int n>
	class characterestic_polynomial
	{
	protected:
		using polynomial_type = std::conditional_t<ring_constraints::field<R>, poly::polynomial<R>, poly::free_algebra<R>>;
	public:
		virtual polynomial_type operator()(const matrix<R, n, n>& A) const = 0;
		auto characterestic(const matrix<R, n, n>& A) const
		{
			return this->operator()(A);
		}
	};
}