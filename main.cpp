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
using K = math_rz::complex;
using E = K;
using E2 = math_rz::linalg::finite_dimensional_vector_space<K, 2>;
using F = K;
using M = math_rz::linalg::square_matrix<K, 2>;
using R_X = math_rz::poly::polynomial<K>;

class mat_exp :public math_rz::analysis::function<E,F>
{
public:
	F operator()(const E& s) const override
	{
		static M rot({ {1,-1},{4,2} });
		return s;
	}

	bool is_zero() const override
	{
		return false;
	}
};

int main()
{
	mat_exp f;
	std::shared_ptr<integrator<K, K>> I_ptr(new trapezoidal_integrator<K, K>(0, 10, 500));
	std::shared_ptr<integrator<K, K>> J_ptr(new trapezoidal_integrator<K, K>(2.+-100.i,2.+ 100.i, 2000));
	laplace_transform L(I_ptr);
	inverse_laplace_transform L_(J_ptr);
	cout << L(f)(1) << L(f)(2) << endl;
	cout << L_(L(f))(5);
	std::shared_ptr<integrator<real_field, real_field>> R_ptr(new trapezoidal_integrator<real_field,real_field >(0, 10, 500));
	std::shared_ptr<integrator<real_field, real_field>> W_ptr(new trapezoidal_integrator<real_field,real_field>(2. , 2., 200));
	cout << R_ptr->integrate(general_function<real_field, real_field>([](const auto& a) {return a; }));
	return false;
}