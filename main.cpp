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
#include "linalg/matrix/ones.h"
#include "linalg/decomposer/polar_decomposition.h"
#include "linalg/decomposer/SVD_decomposition.h"
#include "linalg/matrix/DFT_matrix.h"
using namespace math_rz;
using namespace math_rz::linalg;
using namespace math_rz::analysis;
using K = math_rz::real_field;
using E3 = math_rz::linalg::finite_dimensional_vector_space<K, 3>;
using E2 = math_rz::linalg::finite_dimensional_vector_space<K, 2>;
template<int n>
using E = math_rz::linalg::coordinate_space<K, n>;
using F = K;
using M = math_rz::linalg::matrix<K, 3,5>;
constexpr int dimension = 3;
#include <fstream>
int main()
{
	using namespace std::complex_literals;
	inverter::moore_penrose_pseudo_inverter < K, dimension, dimension > PP;
	special::ones<K, dimension, dimension> A;
	decomposer::polar_decomposition<K, dimension, dimension> PD;
	diagonalisation::QR_algorithm_hermitian<K, dimension> QR;
	special::DFT_matrix<dimension> DFT;
	axes_dilation<E<3>> AD({2,3,4});
	std::cout << DFT << "\n\n";
	//diagonalisation::gram_schmidt_diagonalisation<K, dimension> GSD;
	decomposer::SVD_decomposition<K, dimension, dimension> PDD;
	decomposer::QR_decomposition<K,dimension,dimension> LQ;
	auto [X, Y, Z] = PDD.decompose(A);
	std::cout <<  X <<"\n\n" << Y <<"\n\n" << Z;
}