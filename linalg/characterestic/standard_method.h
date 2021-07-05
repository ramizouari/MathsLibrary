#pragma once
#include "characterestic_polynomial.h"
namespace math_rz::linalg::characterestic
{
	template<typename R, int n>
	class standard_method:public characterestic_polynomial<R,n>
	{
	protected:
		using polynomial_type = std::conditional_t<ring_constraints::field<R>, poly::polynomial<R>, poly::free_algebra<R>>;
	public:
		virtual polynomial_type operator()(const matrix<R, n, n>& A) const
		{
			return A.caracteristic_polynomial();
		}
	};
}