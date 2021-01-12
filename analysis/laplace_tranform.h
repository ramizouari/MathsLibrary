#pragma once
#pragma once
#include "function.h"
#include "complex.h"
#include <numbers>
#include "general_function.h"
#include "integrator/integrator.h"

namespace math_rz
{
	using complex_function_space = function<complex, complex>;
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
}