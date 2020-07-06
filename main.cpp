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
	square_matrix<real_field, 3> P({ {1,1,2},{5,1,3},{2,4,2} });
	coordinate_space<real_field, 3> u({ 1,2,3 });
	cout << rational_extension<polynomial<real_field>>(polynomial<real_field>({1,2}), polynomial<real_field>({4,6}));
	return false;
}