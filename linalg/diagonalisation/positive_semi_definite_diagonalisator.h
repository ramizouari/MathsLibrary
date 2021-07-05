#pragma once
#include "diagonalisator.h"
namespace math_rz::linalg::diagonalisation
{
	template<typename K,int n>
	class positive_semi_definite_diagonalisator:virtual public self_adjoint_diagonalisator<K,n>
	{
		
	};
}