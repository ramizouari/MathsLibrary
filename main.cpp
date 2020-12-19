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
#include <cmath>
#include "analysis/laplace_tranform.h"

using namespace std;
using namespace math_rz;

using K = math_rz::complex;
using E = math_rz::complex;
using F = math_rz::complex;
class sinc_f :public math_rz::function<E,F>
{
	F operator()(const E &s) const override
	{
		if (s.real()<0)
			return 0;
		return 1;
	}

	bool is_zero() const override
	{
		return false;
	}
};


int main()
{
	sinc_f sinc;
	trapezoid_integrator<E,F> I(0,100,2000);
	laplace_transform F(I);
	cout << F(sinc)(1.+1i)<< F(sinc)(5) << endl;
	return false;
}