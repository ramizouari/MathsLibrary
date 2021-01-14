#pragma once
#pragma once
#include "function.h"
#include "complex.h"
#include <numbers>
#include "general_function.h"
#include "integrator/integrator.h"

namespace math_rz::analysis
{
	using complex_function_space = function<complex, complex>;
	/*
	* This class represents the laplace transform of a function
	* It acts with the help of an integrator to find the Laplace Tranform
	* 
	* Constraints on the integrator: the integrator must be a single variable integrator on the positive real time
	* with respect to the standard mesure
	*/
	class laplace_transform : public function<complex_function_space, general_function<complex, complex>>
	{
	public:

		laplace_transform(integrator<complex, complex>* _I_ptr) :I_ptr(_I_ptr) {}
		laplace_transform(std::shared_ptr<integrator<complex, complex>> _I_ptr) :I_ptr(_I_ptr) {}
		general_function<complex, complex> operator()(const complex_function_space& f) const
		{
			auto& J = *I_ptr;
			return general_function<complex, complex>
				(
					[&J, &f](const complex& s)
					{
						return J.integrate
						(
							general_function<complex, complex>([&s, &f](const complex& t)
								{
									return std::exp(-s * t) * f(t);
								})
						);
					}
			);
		}


		/*
		* The laplace transform is not zero on the space of operators:
		*/
		bool is_zero() const override
		{
			return false;
		}
	private:
		std::shared_ptr<integrator<complex, complex>> I_ptr;
	};


	/*
	* This class represents the laplace transform of a function
	* It acts with the help of an integrator to find the Laplace Tranform
	*
	* Constraints on the integrator: the integrator must be a single variable integrator on the line (c-i,c+i) for a well
	* chosen real constant c with respect to the standard mesure
	* The real constant c must be greater then all real component of singularities of the function to be transformed
	*/
	class inverse_laplace_transform : public function<complex_function_space, general_function<complex, complex>>
	{
	public:

		inverse_laplace_transform(integrator<complex, complex>* _I_ptr) :I_ptr(_I_ptr) {}
		inverse_laplace_transform(std::shared_ptr<integrator<complex, complex>> _I_ptr) :I_ptr(_I_ptr) {}
		general_function<complex, complex> operator()(const complex_function_space& f) const
		{
			auto& J = *I_ptr;
			return general_function<complex, complex>
				(
					[&J, &f](const complex& s)
					{
						return J.integrate
						(
							general_function<complex, complex>([&s, &f](const complex& t)
								{
									return std::exp(s * t) * f(t)/complex(0,2*std::numbers::pi);
								})
						);
					}
			);
		}


		/*
		* The inverse laplace transform is not zero on the space of operators:
		*/
		bool is_zero() const override
		{
			return false;
		}
	private:
		std::shared_ptr<integrator<complex, complex>> I_ptr;
	};
}