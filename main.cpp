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

using namespace std;
using namespace math_rz;

using K = math_rz::real_field;
using E = math_rz::Lp_finite_dimensional_space<K,2,2>;
using F = math_rz::real_field;
using M = math_rz::square_matrix<K, 4>;
class mat_exp :public math_rz::function<E,F>
{
public:
	F operator()(const E &s) const override
	{
		return std::pow(s[0],2)*std::pow(s[1],3);
	}

	bool is_zero() const override
	{
		return false;
	}
};


int main()
{
	M H({ {0.4,0.1,0.3,0.2},{0.3,0.4,0.2,0.1},{0.2,0.3,0.1,0.4},{0.1,0.2,0.4,0.3} });
	cout << math_rz::pow(H, 30).caracteristic_polynomial();
	return false;
}