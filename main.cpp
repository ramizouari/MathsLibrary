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
#include <analysis/finite_dimensional_inner_product_space.h>
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
#include "boost/multi_array.hpp"
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

using namespace std;
using namespace math_rz;
using namespace math_rz::linalg;
using namespace math_rz::analysis;
using K = math_rz::real_field;
using E3 = math_rz::linalg::finite_dimensional_vector_space<K, 3>;
using E2 = math_rz::linalg::finite_dimensional_vector_space<K, 2>;
template<int n>
using E = math_rz::linalg::coordinate_space<K, n>;
using F = K;
using M = math_rz::linalg::square_matrix<K, 3>;
using R_X = math_rz::poly::polynomial<K>;

class mat_exp :public math_rz::analysis::function<E<3>,E<3>>
{
public:
	E<3> operator()(const E<3>& s) const override
	{
		static M rot({ {0,-1,0},{1,0,0},{0,0,1} });
		return rot*s;
	}

	bool is_zero() const override
	{
		return false;
	}
};

int main()
{
	uniform_int_generator G(-5, 5, 200);
	constexpr int n = 3000,p=3000,m=3000;
	using matrix_type1 = std::conditional_t<n == p, square_matrix<K, n>, matrix<K, n, p>>;
	using matrix_type2 = std::conditional_t<m == p, square_matrix<K, p>, matrix<K, p, m>>;
	using matrix_type3 = std::conditional_t<n == m, square_matrix<K, n>, matrix<K, n, m>>;
	matrix_type1 A(G.generate_matrix<n,p>());
	matrix_type2 B(G.generate_matrix<p, m>());
	parallel_strassen_multiplicator<K> M(1);
	std::chrono::time_point<std::chrono::system_clock> t1(std::chrono::system_clock::now());
	auto H1=M.multiply<n, p, m>(A, B);
	std::chrono::time_point<std::chrono::system_clock> t2(std::chrono::system_clock::now());
	cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << endl;
	auto H2=(A * B);
	std::chrono::time_point<std::chrono::system_clock> t3(std::chrono::system_clock::now());
	cout << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << endl;
	cout << H2.metric(H1);
	return false;
}