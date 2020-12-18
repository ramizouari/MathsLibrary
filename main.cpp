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
#include <cmath>

using namespace std;

class arctan_fun :public rr_integrable_function
{
	real_field operator()(const real_field &s) const override
	{
		return std::atan(s);
	}

	bool is_zero() const override
	{
		return false;
	}
};


int main()
{
	rectangle_integrator I(0, 1);
	arctan_fun s;
	cout << s.integral(I);
	return false;
}