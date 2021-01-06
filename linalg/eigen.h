#pragma once
#include <random>
#include "square_matrix.h"
#include "analysis/finite_dimensional_inner_product_space.h"
namespace math_rz
{
	template<typename F,int n>
	F largest_eig(const square_matrix<F, n>& M,const real_field &eps=1e-5)
	{
		std::random_device dev;
		std::mt19937 m(dev());
		std::uniform_real_distribution<long double> dist(-1, 1);
		Lp_finite_dimensional_space<F,2, n> u,v,w;
		for (int i = 0; i < n; i++)
			u[i] = dist(m);
		v = M * u;
		if (v.norm()==0)
			return 0;
		v /= v.norm();
		while ((w=(v-u)).norm() > eps)
		{
			u = v;
			v = M * v;
			v /= v.norm();
		}
		int k = 0;
		for(int i=0;i<n;i++)
			if (v[i].norm() > v[k].norm())
				k=i;
		if (v[k].norm() < eps)
			return 0;
		return (M*v)[k]/v[k];
	}
	template<typename F,int n,int m>
	real_field largest_sing(const matrix<F, n,m>& M, const real_field& eps = 1e-5)
	{
		return std::sqrt(static_cast<real_field>(largest_eig(square_matrix<F,n>(M * M.conj_transpose()))));
	}
}