#pragma once
#include "integral_transform.h"
#include <cmath>
#include <numbers>
/*
* This file containts a list of known integral transforms
*/
namespace math_rz::analysis::transforms
{
	class hermite_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
	public:

		hermite_transform(
std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : integral_transform<real_field, real_field, real_field>
			(new general_function<kernel_domain,real_field>([&](const kernel_domain& w)->real_field {return std::hermitel((int)w[1], w[0]); }), _I_ptr) {}
		hermite_transform(integrator<real_field, real_field>* _I_ptr) :
			integral_transform<real_field, real_field, real_field>
			(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field {return std::hermitel((int)w[1], w[0]); }), _I_ptr) {}
	};

	class hartley_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
	public:

		hartley_transform(
			std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : integral_transform<real_field, real_field, real_field>
			(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field 
				{
					real_field p = w[0] * w[1];
					return (std::cos(p)+std::sin(p))/ std::sqrt(2 * std::numbers::pi);;
				}), _I_ptr) {}
		hartley_transform(integrator<real_field, real_field>* _I_ptr) :
			integral_transform<real_field, real_field, real_field>
			(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field 
				{	real_field p = w[0] * w[1];
		return (std::cos(p) + std::sin(p)) / std::sqrt(2*std::numbers::pi);  }), _I_ptr) {}
	};

	class legendre_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
	public:

		legendre_transform(
			std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : integral_transform<real_field, real_field, real_field>
			(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
				{
					return std::legendre((int)w[1], w[0]);
				}), _I_ptr) {}
				legendre_transform(integrator<real_field, real_field>* _I_ptr) :
					integral_transform<real_field, real_field, real_field>
					(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
						{	return std::legendre((int)w[1], w[0]);  }), _I_ptr) {}
	};

	class poisson_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
	public:

		poisson_transform(
			std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : integral_transform<real_field, real_field, real_field>
			(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
				{
					return (1-w[0]*w[0])/(1-2*w[0]*w[1]+w[0]*w[0]);
				}), _I_ptr) {}
				poisson_transform(integrator<real_field, real_field>* _I_ptr) :
					integral_transform<real_field, real_field, real_field>
					(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
						{	
							return (1 - w[0] * w[0]) / (1 - 2 * w[0] * w[1] + w[0] * w[0]);  
						}), _I_ptr) {}
	};


	class weirstrass_transform :public integral_transform<real_field, real_field, real_field>
	{
		using integral_transform<real_field, real_field, real_field>::kernel_domain;
	public:

		weirstrass_transform(
			std::shared_ptr<integrator<real_field, real_field>> _I_ptr) : integral_transform<real_field, real_field, real_field>
			(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
				{
					return std::exp(-std::pow(w[1]-w[0],2)/4)/std::sqrt(4*std::numbers::pi);
				}), _I_ptr) {}
				weirstrass_transform(integrator<real_field, real_field>* _I_ptr) :
					integral_transform<real_field, real_field, real_field>
					(new general_function<kernel_domain, real_field>([&](const kernel_domain& w)->real_field
						{
							return std::exp(-std::pow(w[1] - w[0], 2) / 4) / std::sqrt(4 * std::numbers::pi);
						}), _I_ptr) {}
	};

}