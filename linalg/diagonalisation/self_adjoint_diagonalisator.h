#pragma once
#include "diagonalisator.h"
namespace math_rz::linalg::diagonalisation
{
	template<typename K, int n>
	class self_adjoint_diagonalisator :virtual public diagonalisator<K, n>
	{

	};
}