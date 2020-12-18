#include "integer.h"
#include "guassian_integer.h"
#include "polynomial.h"
#include "ring.h"
#include <iostream>
#include "square_matrix.h"
#include "real_field.h"
#include "finite_dimensional_vector_space.h"
#include "rational_extension.h"
#include "complex.h"

using namespace std;

int main()
{
	square_matrix<integer, 2> M({ {0,1},{1,1} });
	cout << (algebra::pow(M, 40) * finite_dimensional_vector_space<integer, 2>({0,1})).at(0);
	return false;
}