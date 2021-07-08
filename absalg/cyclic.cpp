#include "cyclic.h"
#include <set>
#include <map>
#include <cmath>
math_rz::dynamic_cyclic::dynamic_cyclic(integer _n, integer s):n(_n),w((n._v+s._v)%n)
{
}

using natural = std::uint64_t;
std::set<natural> prime_factors(natural s)
{
	std::set<natural> S;
	natural r = std::ceil(std::sqrt(s));
	for (int i = 2; i <= r;)
		if (s % i == 0)
		{
			s /= i;
			r = std::ceil(std::sqrt(s));
			S.insert(i);
		}
		else i++;
	if(s>1)
		S.insert(s);
	return S;
}

std::map<natural,natural> prime_factors_multiplicity(natural s)
{
	std::map<natural,natural> S;
	natural r = std::ceil(std::sqrt(s));
	for (int i = 2; i <= r;)
		if (s % i == 0)
		{
			s /= i;
			r = std::ceil(std::sqrt(s));
			S[i]++;
		}
		else i++;
	if (s > 1)
		S[s]++;
	return S;
}
natural totient(natural s)
{
	auto P=prime_factors_multiplicity(s);
	natural r = 1;
	for (auto [p, m] : P)
		r *= std::pow(p,m-1)*(p-1);
	return r;
}

using namespace math_rz;
dynamic_cyclic& math_rz::dynamic_cyclic::operator+=(const dynamic_cyclic& b)
{
	if (n != b.n)
		throw std::domain_error("incompatible modulus");
	w += b.w;
	w %= n;
	return *this;
}

dynamic_cyclic& math_rz::dynamic_cyclic::operator-=(const dynamic_cyclic& a)
{

	if (n != a.n)
		throw std::domain_error("incompatible modulus");
	w -= a.w;
	w %= n;
	return *this;
}

dynamic_cyclic& math_rz::dynamic_cyclic::operator*=(const dynamic_cyclic& a)
{
	if (n != a.n)
		throw std::domain_error("incompatible modulus");
	w *= a.w;
	w %= n;
	return *this;
}

dynamic_cyclic& math_rz::dynamic_cyclic::operator+=(int b)
{
	w += b;
	w %= n;
	return *this;
}

dynamic_cyclic& math_rz::dynamic_cyclic::operator-=(int a)
{
	w -=a;
	w %= n;
	return *this;
}

dynamic_cyclic& math_rz::dynamic_cyclic::operator*=(int a)
{
	w *= a;
	w %= n;
	return *this;
}

dynamic_cyclic math_rz::dynamic_cyclic::operator-() const
{
	return dynamic_cyclic(*this,-w);
}

math_rz::dynamic_cyclic::operator integer& ()
{
	return w;
}

math_rz::dynamic_cyclic::operator const integer& () const
{
	return w;
}

math_rz::dynamic_cyclic::operator long long() const
{
	return w;
}

bool math_rz::dynamic_cyclic::is_zero() const
{
	return w==0;
}

bool math_rz::dynamic_cyclic::is_one() const
{
	return w==1;
}

dynamic_cyclic math_rz::dynamic_cyclic::inv() const
{
	return dynamic_cyclic(*this,bezout<integer>(w,n).first);
	//return dynamic_cyclic(1,0);
}

dynamic_cyclic math_rz::dynamic_cyclic::operator/=(const dynamic_cyclic& a)
{
	if (n != a.n)
		throw std::domain_error("incompatible modulus");
	return *this*=a.inv();
}

dynamic_cyclic math_rz::dynamic_cyclic::primitive_root(int m)
{
	auto S = prime_factors_multiplicity(m);
	auto t = totient(m);
	auto P = prime_factors(t);
	for (int k = 2; k < m-1; k++)
	{
		if (!pow(dynamic_cyclic(m, k), t).is_one())
			continue;
		bool primitive = true;
		for (auto p : P)
			if (pow(dynamic_cyclic(m, k), t / p).is_one())
			{
				primitive = false;
				break;
			}
		if (primitive)
			return dynamic_cyclic(m,k);
	}
	return dynamic_cyclic(m, 0);
}

dynamic_cyclic math_rz::dynamic_cyclic::primitive_unity_root(int n, int m)
{
	auto S = prime_factors_multiplicity(n);
	auto P = prime_factors(m);
	for (int k = 2; k < n ; k++)
	{
		if (!pow(dynamic_cyclic(n, k), m).is_one())
			continue;
		bool primitive = true;
		for (auto p : P)
			if (pow(dynamic_cyclic(n, k), m / p).is_one())
			{
				primitive = false;
				break;
			}
		if (primitive)
			return dynamic_cyclic(n, k);
	}
	return dynamic_cyclic(n, 0);
}


dynamic_cyclic math_rz::operator+(const dynamic_cyclic& a, const dynamic_cyclic& b)
{
	auto c = a;
	return c+=b;
}

dynamic_cyclic math_rz::operator-(const dynamic_cyclic& a, const dynamic_cyclic& b)
{
	auto c = a;
	return c -= b;
}

dynamic_cyclic math_rz::operator*(const dynamic_cyclic& a, const dynamic_cyclic& b)
{
	auto c = a;
	return c *= b;
}

dynamic_cyclic math_rz::operator/(const dynamic_cyclic& a, const dynamic_cyclic& b)
{
	auto c = a;
	return c /= b;
}
