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
#include "analysis/integrator/boundary_integrator.h"
#include "analysis/derivator/differential.h"
#include "analysis/integrator/stokes_integrator.h"

using namespace std;
using namespace math_rz;
using namespace math_rz::linalg;
using namespace math_rz::analysis;
using K = math_rz::real_field;
using E3 = math_rz::linalg::finite_dimensional_vector_space<K, 3>;
using E2 = math_rz::linalg::finite_dimensional_vector_space<K, 2>;
template<int n>
using E = math_rz::linalg::coordinate_space<K, n>;
using F = K;
using M = math_rz::linalg::matrix<K, 2>;
using R_X = math_rz::poly::polynomial<K>;

class mat_exp :public math_rz::analysis::function<E<2>,E<2>>
{
public:
	E<2> operator()(const E<2>& s) const override
	{
		static M rot({ {0,-1},{1,0} });
		return rot*s;
	}

	bool is_zero() const override
	{
		return false;
	}
};

int main()
{
	/*
	* Verification of Divergence Theorem
	*/
	mat_exp f;
	/*
	* Integrator along the boundary
	*/
	std::shared_ptr<integrator<real_field, F>> 
		J_ptr_boundary(new multiple_integrator<real_field, F>({0}, 
			{ 2 * std::numbers::pi }, { 2000 }));
	/*
	* Integrator over the manifold
	*/
	std::shared_ptr<integrator<E<2>, F>>
		J_ptr_manifold(new multiple_integrator<E<2>, F>({0, 0 },
			{5, 2 * std::numbers::pi}, {1000, 1000 }));
	/*
	* Derivator over the manifold
	*/
	std::shared_ptr<derivator<E<2>, E<2>>> D_ptr_manifold(new two_way_derivator<E<2>,E<2>>(1e-8));
	/*
	* Derivator along the boundary
	*/
	std::shared_ptr<derivator<real_field, E<2>>> D_ptr_boundary(new two_way_derivator<real_field, E<2>>(1e-8));
	
	general_function<real_field, E<2>> phi_boundary([](const real_field& u) 
		{
			return 5*E<2>({std::cos(u),
				std::sin(u)});
		});

	general_function<E<2>, E<2>> phi_manifold([](const E<2>& u)
		{
			return u[0]*E<2>({ std::cos(u[1]), std::sin(u[1]) });
		});

	differential<E<2>, E<2>> diff(D_ptr_manifold);
	line_integrator<E<2>, E<2>> S(phi_boundary,J_ptr_boundary,D_ptr_boundary);
	disk_integrator<E<2>, F> B(J_ptr_manifold);
	cout << S.integrate((f));
	cout << '\t' << B.integrate(diff.curl(f));
	return false;
}