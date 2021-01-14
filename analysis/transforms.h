#pragma once
#include "integral_transform.h"
#include <cmath>
#include <numbers>
/*
* This file containts a list of some known integral transforms
*/
namespace math_rz::analysis::transforms
{
	/*
	* Hermite Transform
	*/
	class hermite_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
		inline static function<kernel_domain, real_field>* K_ptr = new general_function<kernel_domain, real_field>
			(
				[&](const kernel_domain& w)->real_field 
				{
					return std::hermitel((int)w[1], w[0]); 
				}
			);
	public:

		hermite_transform(std::shared_ptr<integrator<real_field, real_field>> _I_ptr) 
			: integral_transform<real_field, real_field, real_field>(K_ptr, _I_ptr) {}
		hermite_transform(integrator<real_field, real_field>* _I_ptr) :
			integral_transform<real_field, real_field, real_field>(K_ptr, _I_ptr) {}
	};

	/*
	* Hartley Transform
	*/
	class hartley_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
		inline static general_function<kernel_domain, real_field>* K_ptr =
			new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
				{
					real_field p = w[0] * w[1];
					return (std::cos(p) + std::sin(p)) / std::sqrt(2 * std::numbers::pi);;
				});
	public:

		hartley_transform(std::shared_ptr<integrator<real_field, real_field>> _I_ptr) 
			: integral_transform<real_field, real_field, real_field>(K_ptr, _I_ptr) {}
		hartley_transform(integrator<real_field, real_field>* _I_ptr) :
			integral_transform<real_field, real_field, real_field>(K_ptr, _I_ptr) {}
	};

	/*
	* Legendre Transform
	*/
	class legendre_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
		inline static general_function<kernel_domain, real_field>* K_ptr =
			new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
				{
					return std::legendre((long long)w[1], w[0]);
				});
	public:

		legendre_transform(
			std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : integral_transform<real_field, real_field, real_field>
			(K_ptr, _I_ptr) {}
				legendre_transform(integrator<real_field, real_field>* _I_ptr) :
					integral_transform<real_field, real_field, real_field>
					(K_ptr, _I_ptr) {}
	};

	/*
	* Poisson Transform
	*/
	class poisson_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
		inline static general_function<kernel_domain, real_field>* K_ptr =
			new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
				{
					return (1 - w[0] * w[0]) / (1 - 2 * w[0] * w[1] + w[0] * w[0]);
				});
	public:

		poisson_transform(std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : 
			integral_transform<real_field, real_field, real_field>(K_ptr, _I_ptr) {}
		poisson_transform(integrator<real_field, real_field>* _I_ptr) :
			integral_transform<real_field, real_field, real_field>(K_ptr, _I_ptr) {}
	};

	/*
	* Weirtrass Transform
	*/
	class weirstrass_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
		general_function<kernel_domain, real_field>* K_ptr =
			new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
				{
					return std::exp(-std::pow(w[1] - w[0], 2) / 4) / std::sqrt(4 * std::numbers::pi);
				});
	public:

		weirstrass_transform(
			std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : integral_transform<real_field, real_field, real_field>
			(K_ptr, _I_ptr) {}
				weirstrass_transform(integrator<real_field, real_field>* _I_ptr) :
					integral_transform<real_field, real_field, real_field>
					(K_ptr, _I_ptr) {}
	};

	/*
	* Multidimensional Fourier Transform
	*/
	template<linalg::vector_space_constraint::vector_space E>
	class multidimensional_fourier_transform :public integral_transform<E, E, complex>
	{
		using typename integral_transform<E, E, complex>::kernel_domain;
		inline static constexpr int n = E::dimension;
		inline static function<kernel_domain, complex>* K_ptr = new general_function<kernel_domain, complex>([&](const kernel_domain& w)->complex
			{
				E s, t;
				if constexpr (n >= 2) for (int i = 0; i < n; i++)
				{
					t[i] = w[i];
					s[i] = w[i + n];
				}
				else
				{
					t = w[0];
					s = w[1];
				}
				return std::exp(s.dot_product(complex(0, -2 * std::numbers::pi) * t));
			});
	public:

		multidimensional_fourier_transform(
			std::shared_ptr<integrator<E, complex>> _I_ptr) : integral_transform<E, E, complex>
			(K_ptr, _I_ptr) {}
				multidimensional_fourier_transform(integrator<E, complex>* _I_ptr) :
					integral_transform<E, E, complex>(K_ptr, _I_ptr) {}
	};

	/*
	* Multidimensional inverse Fourier Transform
	*/
	template<linalg::vector_space_constraint::vector_space E>
	class multidimensional_inverse_fourier_transform :public integral_transform<E, E, complex>
	{
		using typename integral_transform<E, E, complex>::kernel_domain;
		inline static constexpr int n = E::dimension;
		inline static function<kernel_domain, complex>* K_ptr = new general_function<kernel_domain, complex>([&](const kernel_domain& w)->complex
			{
				E s, t;
				if constexpr (n >= 2) for (int i = 0; i < n; i++)
				{
					t[i] = w[i];
					s[i] = w[i + n];
				}
				else
				{
					t = w[0];
					s = w[1];
				}
				return std::exp(s.dot_product(complex(0, 2 * std::numbers::pi) * t));
			});

	public:

		multidimensional_inverse_fourier_transform(
			std::shared_ptr<integrator<E, complex>> _I_ptr) : integral_transform<E, E, complex>
			(K_ptr, _I_ptr) {}
				multidimensional_inverse_fourier_transform(integrator<E, complex>* _I_ptr) :
					integral_transform<E, E, complex>
					(K_ptr, _I_ptr) {}
	};

	/*
	* Multidimensional Laplace Transform
	*/
	template<linalg::vector_space_constraint::vector_space E>
	class multidimensional_laplace_transform :public integral_transform<E, E, complex>
	{
		using typename integral_transform<E, E, complex>::kernel_domain;
		inline static constexpr int n = E::dimension;
		inline static function<kernel_domain, complex>* K_ptr = new general_function<kernel_domain, complex>([&](const kernel_domain& w)->complex
			{
				E s, t;
				if constexpr (n >= 2) for (int i = 0; i < n; i++)
				{
					t[i] = w[i];
					s[i] = w[i + n];
				}
				else
				{
					t = w[0];
					s = w[1];
				}
				return std::exp(-s.dot_product(t));
			});
	public:

		multidimensional_laplace_transform(
			std::shared_ptr<integrator<E, complex>> _I_ptr) : integral_transform<E, E, complex>
			(K_ptr, _I_ptr) {}
				multidimensional_laplace_transform(integrator<E, complex>* _I_ptr) :
					integral_transform<E, E, complex>(K_ptr, _I_ptr) {}
	};

	/*
	* Multidimensional inverse Laplace Transform
	*/
	template<linalg::vector_space_constraint::vector_space E>
	class multidimensional_inverse_laplace_transform :public integral_transform<E, E, complex>
	{
		using typename integral_transform<E, E, complex>::kernel_domain;
		inline static constexpr int n = E::dimension;
		inline static const complex factor = std::pow<long double>(complex(0, 2 * std::numbers::pi), n);
		inline static general_function<kernel_domain, complex>* K_ptr = new general_function<kernel_domain, complex>([&](const kernel_domain& w)->complex
			{
				E s, t;
				if constexpr (n >= 2) for (int i = 0; i < n; i++)
				{
					t[i] = w[i];
					s[i] = w[i + n];
				}
				else
				{
					t = w[0];
					s = w[1];
				}
				return std::exp(s.dot_product(t)) / factor;
			});
	public:

		multidimensional_inverse_laplace_transform(
			std::shared_ptr<integrator<E, complex>> _I_ptr) : integral_transform<E, E, complex>
			(K_ptr, _I_ptr) {}
				multidimensional_inverse_laplace_transform(integrator<E, complex>* _I_ptr) :
					integral_transform<E, E, complex>
					(K_ptr, _I_ptr) {}
	};

	using fourier_transform = multidimensional_fourier_transform<math_rz::complex>;
	using laplace_transform = multidimensional_laplace_transform<math_rz::complex>;

}