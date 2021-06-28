#include "linalg/matrix.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/inverter/moore_penrose_pseudo_inverter.h"
#include "linalg/decomposer/QR_decomposition.h"
#include "poly/multiplicator/multiplicator.h"
#include "poly/polynomial.h"
#include "poly/structure/inner_product.h"
#include "linalg/matrix/givens_matrix.h"
#include "complex.h"
#include <complex>
#include "linalg/transformation/givens_rotation.h"
#include "linalg/transformation/axes_dilation.h"
#include "linalg/transformation/house_holder_reflection.h"
#include "linalg/transformation/rotation.h"
#include "linalg/diagonalisation/QR_algorithm.h"
#include "linalg/matrix/diagonal.h"
#include "linalg/diagonalisation/gram_schmidt_diagonalisation.h"
using namespace math_rz;
using namespace math_rz::linalg;
using namespace math_rz::analysis;
using K = math_rz::complex;
using E3 = math_rz::linalg::finite_dimensional_vector_space<K, 3>;
using E2 = math_rz::linalg::finite_dimensional_vector_space<K, 2>;
template<int n>
using E = math_rz::linalg::coordinate_space<K, n>;
using F = K;
using M = math_rz::linalg::matrix<K, 3,5>;
#include <fstream>
int main()
{
	using namespace std::complex_literals;
	inverter::moore_penrose_pseudo_inverter < K, 4, 4 > PP;
	matrix<K, 4, 4> A({ {1,1,1,1},{1,2,1,2},{0,0,3,3},{4,1,1,-2} });
	
	diagonalisation::QR_algorithm<K,4> QR;

	std::cout << QR.eigenvalues(A.H()*A) << "\n\n";
	diagonalisation::gram_schmidt_diagonalisation<K, 4> GSD;

	auto [P,D] = GSD.diagonalise(A.H() * A);
	auto [D1,P1] = eigdecomposition(A.H() * A);
	auto [P2,D2] = QR.eigendecomposition(A);
	std::cout << P2 << "\n\n" << D2 << "\n\n" << P2.inv() * A * P2;
}