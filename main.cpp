#include "integer.h"
#include "guassian_integer.h"
#include "polynomial.h"
#include "ring.h"
#include <iostream>
#include "square_matrix.h"
#include "real_field.h"
#include "finite_dimensional_vector_space.h"
#include "rational_extension.h"

using namespace std;

int main()
{
	polynomial<real_field> p({ 0,0,0,0,0,0,1 }), q({ 1,0,2 });
	square_matrix<real_field, 4> P({ {1,2,3,4},{5,6,7,8},{9,10,11,12 },{13,14,15,16}});
	coordinate_space<real_field, 3> u({ 1,2,3 });
	rational_function<real_field> H1(polynomial<real_field>({ 1,2 }), polynomial<real_field>({ 1,1 }));
	rational_function<real_field> H2(polynomial<real_field>({ 1,1 }), polynomial<real_field>({ 1,2 }));
	cout << P.caracteristic_polynomial();
	return false;
}