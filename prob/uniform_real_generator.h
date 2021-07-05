#pragma once
#include "generator.h"
#include "real_field.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/square_matrix.h"
namespace math_rz
{
	class uniform_real_generator :public generator<real_field>
	{
		std::mt19937_64 eng;
		long double a, b;
		long long seed_limit;
		long long counter;
	public:
		uniform_real_generator(long double _a, long double _b, int s_limit) :a(_a), b(_b),
			seed_limit(s_limit), counter(s_limit) {}
		real_field generate()
		{
			static std::random_device rd_dev;
			if (counter >= seed_limit)
			{
				eng.seed(rd_dev());
				counter = 0;
			}
			std::uniform_real_distribution<long double> d(a, b);
			counter++;
			return d(eng);
		}

		template<int n>
		linalg::coordinate_space<real_field, n> generate_vector()
		{
			linalg::coordinate_space<real_field, n> X;
			for (int i = 0; i < n; i++)
				X[i] = generate();
			return X;
		}

		template<int n, int m=n>
		linalg::matrix<real_field, n, m> generate_matrix()
		{
			linalg::matrix<real_field, n, m> X;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					X[i][j] = generate();
			return X;
		}

	};
}