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
using namespace std;
using namespace math_rz;

using K = math_rz::real_field;
using E = math_rz::Lp_finite_dimensional_space<K,2,10>;
using F = math_rz::real_field;
using M = math_rz::square_matrix<real_field, 10>;
using R_X = math_rz::polynomial<K>;
class mat_exp :public math_rz::function<E,F>
{
public:
	F operator()(const E &s) const override
	{
		return std::exp(-std::pow(s.norm(),2)/2);
	}

	bool is_zero() const override
	{
		return false;
	}
};

int main()
{
	L2_finite_dimensional_space<K, 3>A({ 1,2,3 });
	L2_finite_dimensional_space<K, 3>B({ 4,5,15 });
	square_matrix<K, 2> M({ {0,0},{0,0} });
	std::cout << M.rank();
	mat_exp S;
	return false;
}