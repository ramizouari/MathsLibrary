#pragma once
#include "determinant.h"
namespace math_rz::linalg::det
{
	template<ring_constraints::ring R, int n> requires (ring_constraints::field<rational_extension<R>>)
	class guass_method : public determinant<R, n>
	{
	public:
		virtual R det(const matrix<R, n, n>& A) const = 0;
		R operator()(const matrix <R, n, n>& A) const
		{
			using matrix_type_extension = std::conditional_t < field_constraints::field<R>,
				matrix, matrix < rational_extension<R>, n, n >>;
			using ring_extension = std::conditional_t < field_constraints::field<R>,
				R, rational_extension<R>>;
			matrix_type_extension S(*this);
			matrix_type_extension M(S.row_echelon_form());
			ring_extension d(1);
			for (int i = 0; i < n; i++)
				if (M[i][i].is_zero())
					return 0;
				else d *= M[i][i];
			return static_cast<R>(d);
		}
	};
}