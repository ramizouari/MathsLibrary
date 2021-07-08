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
#include "linalg/matrix/centering_matrix.h"
#include "analysis/special/special.h"
#include "linalg/characterestic/interpolation_method.h"
#include "absalg/cyclic.h"
#include "linalg/matrix/vandermonde_matrix.h"
#include "linalg/characterestic/faddev_leverrier.h"
#include "linalg/transformation/fft/cooley_tuckey.h"

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
constexpr int dimension = 3;
constexpr int points=5e6;
#include <fstream>
#include "analysis/integrator/simpson_integrator.h"

int main()
{
	using namespace std::complex_literals;
	matrix<K, dimension> A(special::vandermonde_matrix<K, dimension>({ 1,1,1 }));
	dynamic_cyclic SS(17, 1);
	//fft::cooley_tuckey<8> CT;
	fft::dynamic_cooley_tuckey DCT(dimension);
	characterestic::faddev_leverrier<K, dimension> CP;
	finite_dimensional_vector_space<K, dimension> u,x;
	for (int i = 0; i < dimension; i++)
	{
		u[i] = 1;
		x[i] = i;
	}
	finite_dimensional_vector_space<K, 3> S({ 1,2,3 });
	auto R = poly::fft_interpolation<3>(u);
	cyclic<337> RR;
	finite_dimensional_vector_space<cyclic<337>, 4> T({ 1,2,3,4 });
	fft::dynamic_finite_ring_cooley_tuckey<337> CTR(4);
	fft::fast_convolution FC;
	std::vector<real_field> S1(points);
	for (int i = 0; i < points; i++)
		S1[i] = std::cos(i);
	poly::polynomial<real_field> p(S1);
	std::chrono::time_point t1=std::chrono::system_clock::now();
    decltype(p)::set_multiplicator(new poly::multiplicator::parallel_karatsuba_multiplicator<real_field>());
    auto R1= p * p;
	decltype(p)::set_multiplicator(new poly::multiplicator::fast_multiplicator<real_field>());
	decltype(p)::set_structure(new poly::structure::L2_vect_inner_product<real_field>());
    std::chrono::time_point t2=std::chrono::system_clock::now();
    std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count())/1000 << std::endl;
    auto R2 = p * p;
    std::chrono::time_point t3=std::chrono::system_clock::now();
    std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count())/1000 << std::endl;
    std::cout << (R2 - R1).norm();
	return false;
}