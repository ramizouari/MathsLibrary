#pragma once
#include "generator.h"
#include "absalg/cyclic.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/square_matrix.h"
namespace math_rz
{
	template<int p,bool is_field>
	class uniform_cyclic_generator :public generator<cyclic<p, is_field>>
	{
		std::mt19937_64 eng;
		long long seed_limit;
		long long counter;
	public:
		uniform_cyclic_generator(int s_limit):
			seed_limit(s_limit), counter(s_limit) {}
		cyclic<p, is_field> generate()
		{
			static std::random_device rd_dev;
			if (counter >= seed_limit)
			{
				eng.seed(rd_dev());
				counter = 0;
			}
			std::uniform_int_distribution<long long> d(0, p-1);
			counter++;
			return d(eng);
		}

		template<int n>
		linalg::coordinate_space<cyclic<p, is_field>, n> generate_vector()
		{
			linalg::coordinate_space<cyclic<p, is_field>, n> X;
			for (int i = 0; i < n; i++)
				X[i] = generate();
			return X;
		}

		template<int n, int m=n>
		linalg::matrix<cyclic<p, is_field>, n, m> generate_matrix()
		{
			linalg::matrix<cyclic<p, is_field>, n, m> X;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					X[i][j] = generate();
			return X;
		}
	};
}