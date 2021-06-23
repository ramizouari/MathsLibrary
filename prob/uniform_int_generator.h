#pragma once
#include "generator.h"
#include "integer.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/square_matrix.h"
namespace math_rz
{
	class uniform_int_generator:public generator<integer>
	{
		std::mt19937_64 eng;
		long long a, b;
		long long seed_limit;
		long long counter;
	public:
		uniform_int_generator(long long _a,long long _b,int s_limit):a(_a),b(_b),
			seed_limit(s_limit),counter(s_limit){}
		
		template<typename R=integer>
		R generate()
		{
			static std::random_device rd_dev;
			if (counter >= seed_limit)
			{
				eng.seed(rd_dev());
				counter = 0;
			}
			std::uniform_int_distribution<long long> d(a,b);
			counter++;
			return d(eng);
		}

		integer generate()
		{
			return generate<>();
		}

		template<int n,typename R=integer>
		linalg::coordinate_space<R,n> generate_vector()
		{
			linalg::coordinate_space<R,n> X;
			for (int i = 0; i < n; i++)
				X[i] = generate();
			return X;
		}

		template<int n, int m=n, typename R = integer>
		linalg::matrix<R, n, m> generate_matrix()
		{
			linalg::matrix<R, n, m> X;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					X[i][j] = generate();
			return X;
		}
	};
}