#pragma once
#include "general_function.h"
#include "linalg/finite_dimensional_vector_space.h"
#include "integrator/integrator.h"

namespace math_rz::analysis
{
	template<linalg::vector_space_constraint::vector_space E1,
		linalg::vector_space_constraint::vector_space E2,
		linalg::vector_space_constraint::vector_space F>
	class  integral_transform : public function<function<E1,F>,general_function<E2,F>>
	{
	protected:
		using E = linalg::vector_space_constraint::product_space<E1, E2>;
		//Kernel of the transform
		std::shared_ptr<function<E, F>> K_ptr;
		//Integrator
		std::shared_ptr<integrator<E1, F>> I_ptr;
	public:
		integral_transform(function<E, F>* _K_ptr, integrator<E1, F>*_I_ptr):K_ptr(_K_ptr)
			,I_ptr(_I_ptr){}
		integral_transform(std::shared_ptr<function<E, F>> _K_ptr, integrator<E1, F>* _I_ptr) :K_ptr(_K_ptr)
			, I_ptr(_I_ptr) {}
		integral_transform(function<E, F>* _K_ptr, std::shared_ptr<integrator<E1, F>> _I_ptr) :K_ptr(_K_ptr)
			, I_ptr(_I_ptr) {}
		integral_transform(std::shared_ptr<function<E, F>> _K_ptr, std::shared_ptr<integrator<E1, F>> _I_ptr) 
			:K_ptr(_K_ptr), I_ptr(_I_ptr) {}
		general_function<E2, F> operator()(const function<E1, F>& f) const
		{
			auto& K = *K_ptr;
			return general_function<E2, F>([&](const E2& s)
			{
				return I_ptr->integrate
				(
					general_function<E1, F>([&](const E1& t)
					{
						return f(t) * K(E(t,s));
					})
				);
			});
		}

		bool is_zero() const
		{
			return K_ptr->is_zero();
		}
	};
}