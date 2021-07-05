#pragma once
#include <iostream>
#include <vector>
#include "linalg/finite_dimensional_vector_space.h"

namespace math_rz::multialg 
{
	template<typename K,int ...n>
	class tensor;
	template<typename K>
	class tensor<K>:public K
	{
	public:
		using K::K;

		bool is_zero() const
		{
			return *this == 0;
		}
		static inline constexpr std::vector<int> shape ;
	};
	template<typename K, int n0,int ...n>requires (n0 > 0)
	class tensor<K,n0,n...> : virtual public group
	{

	public:
		tensor() :u(n0)
		{
		}
		tensor(const std::initializer_list<tensor<K,n...>>& M) :
			u(M)
		{
			if (M.size()!=n0)
				throw std::domain_error("Dimensions are not compatible");
		}

	
		tensor(std::initializer_list<tensor<K,n...>>&& M) :
			u([](auto &&S)->auto {
			if (S.size() != n0)
				throw std::domain_error("Dimensions are not compatible");
			return std::move(S); }(std::move(M)))
		{
			
		}

			tensor(const std::vector<tensor<K,n...>>& M) :
			u(M)
		{
			if (M.size()!=n0)
				throw std::domain_error("Dimensions are not compatible");
		}

	
		tensor(std::vector<tensor<K,n...>>&& M) :
			u([](auto &&S)->auto {
			if (S.size() != n0)
				throw std::domain_error("Dimensions are not compatible");
			return std::move(S); }(std::move(M)))
		{
			
		}

		
		inline constexpr static std::vector< int> shape = { n0,n... };
		inline static constexpr int dimension = (n * ...*n0);
		using base_field = K;

		tensor<K, dimension> flatten() const
		{
			if constexpr (sizeof...(n) == 0)
				return *this;
			/*
			* implementation not complete
			*/
			return *this;
		}
		
		bool operator!=(const tensor<K,n0,n...>& M) const
		{
			for (int i = 0; i < n0; i++)
				if (u[i] != M.u[i])
					return true;
			return false;
		}

		bool operator==(const tensor<K,n0,n...>& M) const
		{
			return !(*this != M);
		}
		tensor<K,n0,n...>& operator+=(const tensor<K,n0,n...>& o)
		{
			for (int i = 0; i < n0; i++)
				u[i] += o.u[i];
			return *this;
		}
		tensor<K,n0,n...>& operator-=(const tensor<K,n0,n...>& o)
		{
			for (int i = 0; i < n0; i++)
				u[i] -= o.u[i];
			return *this;
		}
		tensor<K,n0,n...> operator-() const
		{
			tensor<K, n0, n...> S;
			for (int i = 0; i < n0; i++)
				S.u[i] = - u[i];
			return S;
		}

		tensor<K,n0,n...>& operator*=(const K& k)
		{
			for (int i = 0; i < n0; i++)
					u[i] *= k;
			return *this;
		}

		tensor<K,n0,n...>& operator/=(const K& k)
		{
			for (int i = 0; i < n0; i++)
				u[i] /= k;
			return *this;
		}

		auto& operator[](int i)
		{
			return u[i];
		}
		const auto& operator[](int i) const
		{
			return u[i];
		}

		auto& at(int i)
		{
			return u.at(i);
		}
		const auto& at(int i) const
		{
			return u.at(i);
		}

		linalg::finite_dimensional_vector_space<K,n0> 
			operator()(const linalg::finite_dimensional_vector_space<K,n>& u...) const
		{

		}

		bool is_zero() const override
		{
			for (int i = 0; i < n0; i++)
				if (!u[i].is_zero())
					return false;
			return true;
		}

		void foreach(const std::function<void(K&)>& f)
		{
			for (auto& S : u)
				S.foreach(f);
		}
	protected:
		//boost::multi_array<K,2> u;
		//using structure_type = math_rz::multialg::structure::matrix::L22_operator_norm;
		std::vector<tensor<K,n...>> u;

	
	};
	template <typename K, int ...n>
	tensor< K,n...>operator+(
		const tensor< K, n...>& a, const tensor< K, n...>& b)
	{
		auto c(a);
		return c += b;
	}

	template <typename K, int ...n>
	tensor< K, n...>operator-(
		const tensor< K, n...>& a, const tensor< K, n...>& b)
	{
		auto c(a);
		return c -= b;
	}


	template<typename K, int ...n, int p, int ...m>
	tensor<K, n..., m...> operator*(const tensor<K, n...,p>& A, const tensor<K, p, m...>& B)
	{
		tensor<K, n..., m...> C;
		return C;
	}

	template <typename K, int ...n>
	tensor<K, n...> operator*(
		const K& k, const tensor<K,n...>& a)
	{
		auto c(a);
		return c *= k;
	}
	template <typename K, int ...n>
	tensor<K, n...> operator/(
		const K& k, const tensor<K, n...>& a)
	{
		auto c(a);
		return c /= k;
	}

	template <typename K, int n0, int ...n>
	std::ostream& operator<<(std::ostream& H, const tensor<K, n0, n...>& T)
	{
		H << '[';
			for (int i = 0; i < n0; i++)
				H << T[i] << (i < n0-1 ? ", " : "]");
		return H;
	}

	template <typename K, int n0, int ...n>
	std::istream& operator>>(std::istream& H, tensor<K, n0, n...>& T)
	{
		for(int i=0;i<n0;i++)
			H >> T[i];
		return H;
	}
}