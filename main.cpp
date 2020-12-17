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
	square_matrix<::complex, 3> M({ {::complex(1),::complex(2),
		::complex (-1)},{::complex (-2),::complex(0),::complex(1)},
		{::complex(1),::complex (-1),::complex(0)} });
	cout << M.caracteristic_polynomial() << '\n' << (M*M.inv());
	return false;
}