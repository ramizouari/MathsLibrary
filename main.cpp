#include "integer.h"
#include "guassian_integer.h"
#include "poly/polynomial.h"
#include "absalg/ring.h"
#include <iostream>
#include "linalg/square_matrix.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "absalg/rational_extension.h"
#include "complex.h"
#include "analysis/normed_finite_dimensional_space.h"
#include "analysis/finite_dimensional_inner_product_space.h"
#include "analysis/integrable_function.h"
#include "analysis/integrator/integrator.h"
#include "analysis/derivator/derivator.h"
#include "analysis/derivator/default_derivator.h"
#include "analysis/derivator/two_way_derivator.h"
#include "analysis/fourier_transform.h"
#include "analysis/integrator/circular_integrator.h"
#include "analysis/integrator/simpson_integrator.h"
#include <numbers>
#include <cmath>
#include "analysis/laplace_tranform.h"
#include "analysis/integrator/double_integrator.h"
#include "analysis/integrator/triple_integrator.h"
#include "analysis/integrator/spherical_integrator.h"
#include "analysis/integrator/rectangular_integrator.h"
#include "analysis/special/special.h"
#include "poly/roots.h"
#include "absalg/cyclic.h"
#include "absalg/ring_extension.h"
#include "linalg/eigen.h"
#include "analysis/integrator/multiple_integrator.h"
#include "poly/interpolation.h"
#include "prob/uniform_int_generator.h"
#include "prob/uniform_cyclic_generator.h"
#include "prob/uniform_real_generator.h"
#include "prob/uniform_complex_generator.h"
#include "poly/multiplicator/multiplicator.h"
#include "poly/structure/inner_product.h"
#include "analysis/integrator/disk_integrator.h"
#include "analysis/integrator/ball_integrator.h"
#include "analysis/integrator/line_integrator.h"
#include "analysis/integrator/surface_integrator.h"
#include "analysis/integrator/substitution_integrator.h"
#include "linalg/structure/matrix/inner_product.h"
#include "analysis/structure/function/inner_product.h"
#include "analysis/integral_transform.h"
#include "analysis/transforms.h"
#include "analysis/integrator/boundary_integrator.h"
#include "analysis/derivator/differential.h"
#include "analysis/integrator/stokes_integrator.h"
#include "linalg/multiplicator/multiplicator.h"
#include <chrono>
#include "analysis/minimser/gradient_descent.h"
#include "analysis/minimser/fixed_rate_gradient_descent.h"
#include "analysis/minimser/barzilai_borwein_gradient_descent.h"
#include "analysis/solver/fixed_point_iteration.h"
#include "analysis/solver/newton_raphson.h"
#include "analysis/solver/secant_method.h"
#include "linalg/decomposer/QR_decomposition.h"
#include "linalg/inverter/moore_penrose_pseudo_inverter.h"

using namespace std;
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
using R_X = math_rz::poly::polynomial<K>;

#include <fstream>
int main()
{
	uniform_complex_generator G(0,20,0,20,500);
	matrix<K, 300, 200> A=G.generate_matrix<300,200>();
	E3 b({ 1,1,1 });
	square_matrix<K, 3> I(1);
	std::ofstream file("out.txt");
	decomposer::cholesky<K,300> Ch;
	decomposer::QR_decomposition<K,300,200> QR;
	inverter::moore_penrose_pseudo_inverter<K,300,200> P_INV;
	auto [Q, R] = QR.decompose(A);
	file << P_INV.pinv(A);
	file << "\n\n";
	file << (pow(square_matrix<K,300>(A*P_INV.pinv(A)),2)-A*P_INV.pinv(A)).norm() ;
}