#include <numbers>
#include <iostream>
#include <cmath>
#include "integer.h"
#include "linalg/square_matrix.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "analysis/integrator/integrator.h"
#include "analysis/derivator/default_derivator.h"
#include "analysis/derivator/two_way_derivator.h"
#include "analysis/integrator/multiple_integrator.h"
#include "analysis/integrator/line_integrator.h"
#include "analysis/integrator/surface_integrator.h"
#include "analysis/derivator/differential.h"

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
using M = math_rz::linalg::square_matrix<K, 3>;
using R_X = math_rz::poly::polynomial<K>;

class mat_exp :public math_rz::analysis::function<E<3>, E<3>>
{
public:
	E<3> operator()(const E<3>& s) const override
	{
		static M rot({ {0,-1,0},{1,0,0},{0,0,1} });
		return rot * s;
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
		J_ptr_boundary(new multiple_integrator<real_field, F>({ 0 },
			{ 2 * std::numbers::pi }, { 2000 }));
	/*
	* Integrator over the manifold
	*/
	std::shared_ptr<integrator<E<2>, F>>
		J_ptr_manifold(new multiple_integrator<E<2>, F>({ 0, 0 },
			{ 5, 2 * std::numbers::pi }, { 1000, 1000 }));
	/*
	* Derivator over the manifold
	*/
	std::shared_ptr<derivator<E<2>, E<3>>> D_ptr_manifold(new two_way_derivator<E<2>, E<3>>(1e-8));
	/*
	* Derivator used by differential
	*/
	std::shared_ptr<derivator<E<3>, E<3>>> D_ptr_differential(new two_way_derivator<E<3>, E<3>>(1e-8));

	/*
	* Derivator along the boundary
	*/
	std::shared_ptr<derivator<real_field, E<3>>> D_ptr_boundary(new two_way_derivator<real_field, E<3>>(1e-8));

	general_function<real_field, E<3>> phi_boundary([](const real_field& u)
		{
			return 5 * E<3>({ std::cos(u),
				std::sin(u),0 });
		});

	general_function<E<2>, E<3>> phi_manifold([](const E<2>& u)
		{
			return u[0] * E<3>({ std::cos(u[1]), std::sin(u[1]),0 });
		});

	/*
	* Creating a differential operator acting on the space of vector functions (form R3 to itself)
	*/
	differential<E<3>, E<3>> diff(D_ptr_differential);
	/*
	* Creating a line integrator going along the circle of radius 5 centered at the origin
	* located on the xy axis
	*/
	line_integrator<E<3>, E<3>> S(phi_boundary, J_ptr_boundary, D_ptr_boundary);
	/*
	* Creating a surface integrator going over the disk of radius 5 centered at the origin located on the
	* xy axis
	*/
	surface_integrator<E<3>, E<3>> B(phi_manifold, J_ptr_manifold, D_ptr_manifold);
	cout << S.integrate((f));
	cout << '\t' << B.integrate(diff.curl(f));
	return false;
}