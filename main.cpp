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
#include "multialg/tensor.h"
#include "analysis/minimser/lagrangian_minimiser.h"
#include "analysis/minimser/constrained_barzilai_borwein_gradient_descent.h"
#include "analysis/minimser/karush_kuhn_tucker_minimiser.h"
#include "linalg/diagonalisation/QR_algorithm.h"

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
#include "linalg/matrix/house_holder_matrix.h"
#include "linalg/matrix/circulant_matrix.h"
#include <fstream>
#include "analysis/integrator/simpson_integrator.h"

int main()
{
	using namespace std::complex_literals;
	matrix<K, dimension> A(special::vandermonde_matrix<K, dimension>({ 1,1,1 }));
	dynamic_cyclic SS(17, 1);
	fft::cooley_tuckey<8> CT;
	fft::dynamic_cooley_tuckey DCT(dimension);
	characterestic::faddev_leverrier<K, dimension> CP;
	finite_dimensional_vector_space<K, dimension> u,x;
	for (int i = 0; i < dimension; i++)
	{
		u[i] = 1;
		x[i] = i;
	}
	finite_dimensional_vector_space<K, 3> S({ 1,2,3 });
	auto R = poly::fft_interpolation(u);
	cyclic<337> RR;
	finite_dimensional_vector_space<cyclic<337>, 4> T({ 1,2,3,4 });
	fft::dynamic_finite_ring_cooley_tuckey<337> CTR(4);
	fft::fast_convolution FC;
	std::vector<real_field> S1(5e5);
	for (int i = 0; i < 5e5; i++)
		S1[i] = std::cos(i);
	poly::polynomial<real_field> p(S1);
	auto R1= p * p;
	decltype(p)::set_multiplicator(new poly::multiplicator::fast_multiplicator<real_field>());
	decltype(p)::set_structure(new poly::structure::L2_vect_inner_product<real_field>());
	auto R2 = p * p;
	std::cout << (R2 - R1).norm();
	return false;
}