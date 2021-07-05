#pragma once
#include "linalg/matrix.h"
namespace math_rz::linalg::det
{
	template<typename R,int n>
	class determinant
	{
	public:
		virtual R det(const matrix<R, n, n>& A) const = 0;
		R operator()(const matrix <R, n, n>& A) const
		{
			return det(A);
		}
	};
}