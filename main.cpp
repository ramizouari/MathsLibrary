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

using namespace std;
using namespace math_rz;
using K = math_rz::real_field;
using E = math_rz::L2_finite_dimensional_space<real_field,2>;
using F = E;
using M = math_rz::square_matrix<real_field, 10>;
using R_X = math_rz::polynomial<K>;
class mat_exp :public math_rz::function<E,F>
{
public:
	F operator()(const E &s) const override
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
	trapezoidal_integrator<real_field, real_field>*I(
		new trapezoidal_integrator<real_field, real_field>(0,2*std::numbers::pi,100));
	auto f = general_function<real_field, E>([](const real_field& a)->E
		{return E({ std::cos(a),std::sin(a) }); });
	default_derivator<real_field, real_field::dimension, E::dimension, real_field, E>*
		D(new default_derivator<real_field, real_field::dimension, E::dimension, real_field, E>(0,1e-5));
	line_integrator<E, F> L(f,I,D);
	mat_exp w;
	cout << L.integrate(w);
	return false;
}