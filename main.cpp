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
#include "analysis/integrator.h"
#include "analysis/derivator.h"
#include <cmath>

using namespace std;
using namespace math_rz;

using E = finite_dimensional::hilbert_space<3>;
class arctan_fun :public function<E,E>
{
	finite_dimensional::hilbert_space <3> operator()(const E &s) const override
	{
		return E({s.at(0),s.at(1),s.at(2)});
	}

	bool is_zero() const override
	{
		return false;
	}
};


int main()
{
	rectangle_integrator<finite_dimensional::hilbert_space<3>> I(0, 1);
	trapezoid_integrator<finite_dimensional::hilbert_space<3>> J(0, 1);
	arctan_fun s;
	derivator<math_rz::complex, 3, 3> D(math_rz::finite_dimensional::hilbert_space<3>({0,0,0}), .1i);
	cout << D.jacobian(s);
	return false;
}