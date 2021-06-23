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
	inverter::moore_penrose_pseudo_inverter<K,3,3> PP;
	matrix<K, 3> A({ {1,2,3},{1,2,3}, {1,2,3} });
	special::givens_matrix<K,3> B(0,1,1.+5.i,1);
	finite_dimensional_vector_space<K, 3> u({ 1,2,5 });
	givens_rotation<decltype(u)> GR(0,1, 1. + 5.i, 1);
	axes_dilation<decltype(u)> AD({ 2,1,3 });
	house_holder_reflection<decltype(u)> HHR({ 0,0,1 });
	poly::polynomial<K>::set_structure(new poly::structure::L2_vect_inner_product<K>);
	std::cout << B*u << "\n\n" << HHR(u);
}