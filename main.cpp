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

using namespace std;
using namespace math_rz;
using namespace math_rz::linalg;
using namespace math_rz::analysis;
using K = math_rz::real_field;
using E = math_rz::linalg::finite_dimensional_vector_space<real_field,3>;
using F = real_field;
using M = math_rz::linalg::square_matrix<math_rz::complex, 3>;
using R_X = math_rz::poly::polynomial<K>;

class mat_exp :public math_rz::analysis::function<E,F>
{
public:
	K operator()(const E& s) const override
	{
		static M rot({ {1,-1,1},{2,6,3},{5,2,4} });
		return s.norm();
	}

	bool is_zero() const override
	{
		return false;
	}
};


class function_class :public math_rz::analysis::function<finite_dimensional_vector_space<K,2>, K>
{
public:
	K operator()(const finite_dimensional_vector_space<K,2>& s) const override
	{
		return s[0];
	}

	bool is_zero() const override
	{
		return false;
	}
};

class id :public math_rz::analysis::function<K, K>
{
public:
	K operator()(const K& s) const override
	{
		return s;
	}

	bool is_zero() const override
	{
		return false;
	}
};

int main()
{
	square_matrix<math_rz::complex, 3> H({ {59. + 1.i,5,1},{62. + 5.i,6,1.i},{2,8,1} }), S = H.T();
	uniform_cyclic_generator<2, true> G(200);
	square_matrix<cyclic_field<2>, 512> W(G.generate_matrix<512>());
	using namespace finite_dimensional;
	using V = finite_dimensional_vector_space<real_field, 2>;
	std::shared_ptr<rectangular_integrator<real_field, real_field>> I ( new
		rectangular_integrator<real_field, real_field>(-10, 10, 100));
	mat_exp* N = new mat_exp;
	integral_transform<K,V , real_field> T(N, I);
	function_class f;
	id x;
	transforms::hartley_transform hart(I);
	cout << hart(hart(x))(1) << hart(hart(x))(2);
	//cout << T(x)(1) << T(f)(V({2, 2}));
	return false;
}