#include "integer.h"
#include "guassian_integer.h"
#include "polynomial.h"
#include "ring.h"
#include <iostream>
#include "square_matrix.h"
#include "real_field.h"
#include "finite_dimensional_vector_space.h"

using namespace std;

int main()
{
	square_matrix<real_field, 3> P({ {1,1,2},{5,1,3},{2,4,2} });
	coordinate_space<integer, 2> u({ 1,2 });
	cout << P.det() << 6*u;
	return false;
}