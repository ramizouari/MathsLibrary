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

using namespace std;
using namespace math_rz;
using K = math_rz::real_field;
using E = math_rz::L2_finite_dimensional_space<real_field,3>;
using F = E;
using M = math_rz::square_matrix<real_field, 3>;
using R_X = math_rz::polynomial<K>;
class mat_exp :public math_rz::function<E,F>
{
public:
	F operator()(const E &s) const override
	{
		static M rot({ {1,-1,1},{2,6,3},{5,2,4} });
		return (rot * E({ s[0] * s[0],s[1] * s[1],s[2] * s[2]}));
	}

	bool is_zero() const override
	{
		return false;
	}
};

int main()
{
	mat_exp f;
	auto xi = general_function<E, E>
		([](const E& u)->E
			{
				return E({ u[0] * std::cos(u[1])*std::sin(u[2]),
					u[0] * std::sin(u[1])*std::sin(u[2]),
					u[0]*std::cos(u[2]) });
			});
	std::shared_ptr< trapezoidal_triple_integrator<F>> I_ptr
	(new trapezoidal_triple_integrator<F>(0, 2, 0, 2 * std::numbers::pi, 0, std::numbers::pi, 50, 50, 50));
	substitution_integrator<E,F> J(xi, I_ptr,
		new two_way_derivator<E, E>(1e-8));
	ball_integrator<F> B(I_ptr);
	std::cout << J.integrate(f) << endl << B.integrate(f);
	return false;
}