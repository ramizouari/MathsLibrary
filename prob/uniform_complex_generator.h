#pragma once
#include "generator.h"
#include "complex.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/square_matrix.h"
namespace math_rz
{
	class uniform_complex_generator :public generator<complex>
	{
		std::mt19937_64 eng;
		long double a, b,c,d;
		long long seed_limit;
		long long counter;
	public:
		uniform_complex_generator(long double _a, long double _b,long double _c,long double _d, int s_limit) 
			:a(_a), b(_b),c(_c),d(_d),
			seed_limit(s_limit), counter(s_limit) {}
		complex generate()
		{
			static std::random_device rd_dev;
			std::uniform_real_distribution<long double> d1(a, b);
			std::uniform_real_distribution<long double> d2(c, d);
			if (counter >= seed_limit)
			{
				eng.seed(rd_dev());
				counter = 0;
			}
			counter++;
			return complex(d1(eng),d2(eng));
		}

		template<int n>
		linalg::coordinate_space<complex, n> generate_vector()
		{
			linalg::coordinate_space<complex, n> X;
			for (int i = 0; i < n; i++)
				X[i] = generate();
		}

		template<int n, int m>
		linalg::matrix<complex, n, m> generate_matrix()
		{
			linalg::matrix<complex, n, m> X;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					X[i][j] = generate();
			return X;
		}

		template<int n>
		linalg::square_matrix<complex, n> generate_matrix()
		{
			return generate_matrix<n, n>();
		}
	};
}