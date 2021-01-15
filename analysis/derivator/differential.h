#pragma once
#include "derivator.h"
#include "analysis/general_function.h"
#include "linalg/square_matrix.h"
namespace math_rz::analysis
{
	/*
	* This alias gives the type of the differential:
	* Either a matrix or a square matrix
	*/
	template<linalg::vector_space_constraint::normed_vector_space E,
		linalg::vector_space_constraint::normed_vector_space F>
		requires linalg::vector_space_constraint::vector_space_over_same_base_field < E,F>
		using differential_type = std::conditional_t<E::dimension==F::dimension,
		linalg::square_matrix<typename E::base_field,E::dimension>,
		linalg::matrix<typename E::base_field, F::dimension, E::dimension>>;

	/*
	* This class 
	*/
	template<linalg::vector_space_constraint::normed_vector_space E,
		linalg::vector_space_constraint::normed_vector_space F> 
		requires linalg::vector_space_constraint::vector_space_over_same_base_field < E,F>
		class differential :public function<function<E, F>, general_function<E,differential_type<E,F>>>
	{
		std::shared_ptr<derivator<E, F>> D_ptr;
		using K = typename E::base_field;
		using curl_type = std::conditional_t<E::dimension == 3, E, K>;
	public:
		differential(std::shared_ptr<derivator<E,F>> _D_ptr):D_ptr(_D_ptr){}
		differential(derivator<E,F> *_D_ptr):D_ptr(_D_ptr){}

		general_function<E, differential_type<E, F>> operator()(const function<E, F>& f) const override
		{
			return general_function<E,differential_type<E,F>>([&](const E& u)->differential_type<E, F>
			{
				return D_ptr->jacobian(f, u);
			});
		}

		general_function<E, differential_type<E, F>> jacobian(const function<E, F>& f) const
		{
			return this->operator()(f);
		}

		general_function<E,K> jacobian_det(const function<E, F>& f) const
			requires (E::dimension == F::dimension)
		{
			return general_function<E, K>([&](const E& u)->E
			{
				return D_ptr->jacobian(f, u).det();
			});
		}

		general_function<E,E> gradient(const function<E, F>& f) requires (F::dimension == 1)
		{
			return general_function<E,E>([&](const E& u)->E
			{
				return D_ptr->gradient(f, u);
			});
		}

		general_function<E, K> divergence(const function<E, F>& f)const  requires (F::dimension == E::dimension)
		{
			return general_function<E,K>([&](const E& u)->K
			{
				return D_ptr->divergence(f, u);
			});
		}

		general_function<E, curl_type> curl(const function<E, F>& f) const requires (F::dimension == E::dimension) && (E::dimension == 3 || E::dimension==2)
		{
			return general_function<E, curl_type>([&](const E& u)->curl_type
			{
				return D_ptr->curl(f, u);
			});
		}

		/*
		* The differential operator is not zero on the space of smooth functions
		*/
		bool is_zero() const
		{
			return false;
		}
	};
}