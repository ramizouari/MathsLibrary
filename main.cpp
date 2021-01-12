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

using namespace std;
using namespace math_rz;
using K = math_rz::real_field;
using E = math_rz::complex;
using F = E;
using M = math_rz::square_matrix<math_rz::complex, 3>;
using R_X = math_rz::polynomial<K>;

class mat_exp :public math_rz::function<E,F>
{
public:
	F operator()(const E &s) const override
	{
		static M rot({ {1,-1,1},{2,6,3},{5,2,4} });
		if (s.real().abs() < 1)
			return 1;
		return 0;
	}

	bool is_zero() const override
	{
		return false;
	}
};

class mat_exp2 :public math_rz::function<E, F>
{
public:
	F operator()(const E& s) const override
	{
		static M rot({ {1,-1,1},{2,6,3},{5,2,4} });
		return std::exp(-math_rz::complex(std::numbers::pi)*s*s);
	}

	bool is_zero() const override
	{
		return false;
	}
};

int main()
{
	mat_exp f;
	M P({ {1,1,0},{1,-1,1},{2,2,2} });
	cout << P << endl;
	using namespace math_rz::linalg::structure::matrix;
	L22_operator_norm<math_rz::complex, 3, 3> op;
	L2_vect_inner_product<math_rz::complex, 3, 3> op2;
	constexpr auto inf = std::numeric_limits<long double>::infinity();
	Lpq_vector_norm<math_rz::complex, 3, 3> op3(inf, inf),op4(1, inf),op5(inf,1);
	L1q_operator_norm<math_rz::complex,3,3> op6(inf),op7(2);
	M::set_structure(new math_rz::linalg::structure::matrix::L2_vect_inner_product<math_rz::complex,3,3>);
	cout << op.norm(P) << '\t' << op2.norm(P)<<'\t' << op3.norm(P) << '\t' << op4.norm(P) << '\t'
		<< op5.norm(P) << '\t' << op6.norm(P) << '\t' <<  op7.norm(P) << '\t' << P.norm();
	cout << endl << P.inner_product(P)<<endl;
	auto pos_def_matrix = square_matrix<math_rz::complex, 3>({ {2,1,0},{1,1,0},{0,0,1} });
	finite_dimensional_vector_space<math_rz::complex, 3>::set_structure
	(new math_rz::linalg::structure::vector::L2_induced_vect_inner_product<math_rz::complex,3>
		(pos_def_matrix));
	finite_dimensional_vector_space<math_rz::complex, 3>u ({ 1,-1,0 });
	cout << u.norm();
	cout << endl << math_rz::largest_sing(pos_def_matrix)<<endl;
	std::shared_ptr<integrator<math_rz::complex, math_rz::complex>> I
	(new trapezoidal_integrator<math_rz::complex, math_rz::complex>(-10,10,1000)),
		H(new rectangular_integrator<math_rz::complex, math_rz::complex>(-1, 1, 1000));
	laplace_transform L(I);
	inverse_fourier_transform F_(I);
	fourier_transform F(I);
	mat_exp2 g;
	using function_structure = analysis::structure::function::L2_function_inner_product<math_rz::complex, math_rz::complex>;
	function_structure* op_L2(new function_structure(H));
	cout << F(g)(1) <<'\t' << g(1);
	cout << F(g)(0) << '\t' << g(0);
	complex_function_space::set_structure(op_L2);
	cout << endl << (F(g)*F(g)).inner_product(g*g);
	return false;
}