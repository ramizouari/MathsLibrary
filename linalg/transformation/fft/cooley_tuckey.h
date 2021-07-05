#pragma once
#include <cmath>
#include <numeric>
#include "linalg/transformation/linear_transformation.h"
#include <numbers>
#include <type_traits>
#include "linalg/transformation/axes_dilation.h"
#include "absalg/cyclic.h"
namespace math_rz::linalg::fft
{
	template<unsigned int n,bool inverse=false> requires (std::has_single_bit(n)==1)
	class cooley_tuckey:public endomorphism<finite_dimensional_vector_space<complex,n>>
	{
		cooley_tuckey<n / 2,inverse> CT;
		using E = finite_dimensional_vector_space<complex, n>;
		inline static complex w = complex(std::polar<long double>(1,(-2*static_cast<int>(inverse)+1)*2*std::numbers::pi/n));
		axes_dilation<E> AD;
	public:
		cooley_tuckey():AD
			(
				[]()
				{
					E W;
					W[0] = 1;
					for (int i = 1; i < n; i++)
						W[i] = W[i - 1] * w;
					return W;
				}
				()
			)
		{

		}
		E& apply(E& u) const override
		{
			return u=this->operator()(u);
		}

		E operator()(const E& a) const
		{
			using F = finite_dimensional_vector_space<complex, n / 2>;
			F X[2];
			for (int i = 0; i < n; i++)
				X[i % 2][i / 2] = a[i];
			F Y[2] = { CT(X[0]),CT(X[1]) };

			E Z1(Y[0], Y[0]), Z2(Y[1], Y[1]);
			return Z1 + AD.apply(Z2);
		}

		E scaled(const E& a) const
		{
			auto y = this->operator()(a);
			auto N = std::sqrt(n);
			y.foreach([N](auto& x) {x /= N; });
			return y;
		}

		bool is_zero() const
		{
			return false;
		}
	};

	template<bool inverse>
	class cooley_tuckey<1,inverse> :public endomorphism<finite_dimensional_vector_space<complex, 1>>
	{
	public:
		using E = finite_dimensional_vector_space<K, 1>;
		E& apply(E& u) const override
		{
			return u;
		}


		E scaled(const E& u) const
		{
			return u;
		}


		bool is_zero() const
		{
			return false;
		}
	};

	template<int n>
	using reverse_cooley_tuckey = cooley_tuckey<n, true>;

	unsigned long long smallest_prime_divisor(unsigned long long n)
	{
		unsigned long long s = std::ceil(std::sqrt(n));
		for (int p = 2; p <= s; p++)
			if (n % p == 0)
				return p;
		return n;
	}

	template<bool inverse = false>
	class dynamic_cooley_tuckey
	{
		int n;
		complex w;
		using E = std::vector<complex>;
	public:
		dynamic_cooley_tuckey(int _n) :n(_n), w(std::polar<long double>(1,(-2*static_cast<int>(inverse)+1)*2*std::numbers::pi/n))
		{

		}

		std::vector<complex> operator()(const std::vector<complex>& a) const
		{
			if (n == 1)
				return a;
			auto p = smallest_prime_divisor(n);
			auto q = n / p;
			dynamic_cooley_tuckey<inverse> CT(q);
			std::vector<E> X(p,E(q));
			for (int i = 0; i < n; i++)
				X[i % p][i / p] = a[i];
			std::vector<E> Y(p);
			for (int i = 0; i < p; i++)
				Y[i] = CT(X[i]);
			std::vector<complex> W(n);
			W[0] = 1;
			for (int i = 1; i < n; i++)
				W[i] = w * W[i - 1];
			E Z(n);
			for (int i = 0; i < n; i++)
				for(int m=0;m<p;m++)
					Z[i] += Y[m][i%q]*W[(i*m)%n];
			return Z;
		}

	};

	template<unsigned int m,bool is_prime=false>
	class dynamic_finite_ring_cooley_tuckey
	{
		int n;
		using F = cyclic<m,is_prime>;
		F w;
		using E = std::vector<F>;
	public:
		dynamic_finite_ring_cooley_tuckey(int _n) :n(_n), w(F::primitive_unity_root(n))
		{

		}

		std::vector<F> operator()(const std::vector<F>& a) const
		{
			if (n == 1)
				return a;
			auto p = smallest_prime_divisor(n);
			auto q = n / p;
			dynamic_finite_ring_cooley_tuckey<m,is_prime> CT(q);
			std::vector<E> X(p, E(q));
			for (int i = 0; i < n; i++)
				X[i % p][i / p] = a[i];
			std::vector<E> Y(p);
			for (int i = 0; i < p; i++)
				Y[i] = CT(X[i]);
			std::vector<F> W(n);
			W[0] = 1;
			for (int i = 1; i < n; i++)
				W[i] = w * W[i - 1];
			E Z(n);
			for (int i = 0; i < n; i++)
				for (int r = 0; r < p; r++)
					Z[i] += Y[r][i % q] * W[(i * r) % n];
			return Z;
		}
	};
}