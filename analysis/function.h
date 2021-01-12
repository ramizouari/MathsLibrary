#pragma once
#include "linalg/vector_space.h"
#include "absalg/ring.h"
#include "structure/function/inner_product.h"

namespace math_rz {
	template<typename A, typename B>
	class abstract_function
	{
	public:
		using domain = A;
		using codomain = B;
		virtual B operator()(const A& a) const = 0;

	};


	template<typename A, typename B>
	class summable_function : public abstract_function<A, B>, public vector_space<B>
	{
	public:

	};
	template<typename A, typename B>
	class function : public summable_function<A, B>
	{
	protected:
		using structure_type = analysis::structure::function::metric_topology<A, B>;
		inline static std::unique_ptr<structure_type> structure_ptr;
	public:
		static void set_structure(structure_type* S)
		{
			structure_ptr.reset(S);
		}
		static const structure_type& get_structure()
		{
			return (*structure_ptr);
		}

		real_field metric(const function& p)
		{
			return structure_ptr->metric(*this, p);
		}

		real_field distance(const function& p)
		{
			return metric(p);
		}
	};


	template<typename A, math_rz::vector_space_constraint::vector_space B>
	class function<A,B> : public summable_function<A, B>
	{
	protected:
		using structure_type = analysis::structure::function::metric_topology<A, B>;
		using K = typename analysis::structure::function::inner_product_topology<A, B>::K;
		inline static std::unique_ptr<structure_type> structure_ptr;
	public:
		static void set_structure(structure_type* S)
		{
			structure_ptr.reset(S);
		}
		static const structure_type& get_structure()
		{
			return (*structure_ptr);
		}

		real_field metric(const function& p)
		{
			return structure_ptr->metric(*this, p);
		}

		real_field distance(const function& p)
		{
			return metric(p);
		}

		real_field norm() const
		{
			return dynamic_cast<math_rz::analysis::structure::function::norm_topology<A, B>*>
				(structure_ptr.get())->norm(*this);
		}
		K inner_product(const function& q) const
		{
			return dynamic_cast<math_rz::analysis::structure::function::inner_product_topology<A,B>*>
				(structure_ptr.get())->inner_product(*this, q);
		}
	};

	template<typename A,typename B>
	class function_sum:public function<A, B>
	{
		const function<A, B> &f, &g;
	public:
		function_sum(const function<A,B>&_f,const function<A,B>&_g):f(_f),g(_g){}
		virtual B operator()(const A& a) const override
		{
			return f(a) + g(a);
		}
		bool is_zero() const
		{
			return false;
		}
	};

	template<typename A, typename B>
	class function_difference:public function<A,B>
	{
		const function<A, B>& f, & g;
	public:
		function_difference(const function<A, B>& _f, const function<A, B>& _g) :f(_f), g(_g) {}
		virtual B operator()(const A& a) const override
		{
			return f(a) - g(a);
		}
		bool is_zero() const
		{
			return false;
		}
	};


	template<typename A, typename B>
	class function_product:public function<A,B>
	{
		const function<A, B>& f, & g;
	public:
		function_product(const function<A, B>& _f, const function<A, B>& _g) :f(_f), g(_g) {}
		virtual B operator()(const A& a) const
		{
			return f(a) * g(a);
		}
		bool is_zero() const
		{
			return false;
		}
	};

	template<typename A, typename B>
	class function_division:public function<A,B>
	{
		const function<A, B>& f, & g;
	public:
		function_division(const function<A, B>& _f, const function<A, B>& _g) :f(_f), g(_g) {}
		virtual B operator()(const A& a) const
		{
			return f(a) / g(a);
		}
		bool is_zero() const
		{
			return false;
		}
	};

	template<typename A,typename B>
	function_sum<A, B> operator+(const function<A, B>& a, const function<A, B>& b)
	{
		return function_sum<A,B>(a, b);
	}

	template<typename A, typename B>
	function_difference<A, B> operator-(const function<A, B>& a,const  function<A, B>& b)
	{
		return function_difference<A,B>(a, b);
	}

	template<typename A, typename B>
	function_product<A, B> operator*(const function<A, B>& a, const function<A, B>& b)
	{
		return function_product<A, B>(a, b);
	}

	template<typename A, typename B>
	function_division<A, B> operator/(const function<A, B>& a, const  function<A, B>& b)
	{
		return function_division<A, B>(a, b);
	}
}