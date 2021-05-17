#include "real_field.h"

using namespace math_rz;
//real_field::real_field(double m) :_v(m)
//{
//}
real_field real_field::abs() const
{
	return std::abs(_v);
}

real_field math_rz::real_field::norm() const
{
	return abs();
}

real_field math_rz::real_field::conj() const
{
	return _v;
}

real_field math_rz::real_field::inner_product(const real_field& a)const
{
	return (*this) * a;
}

real_field math_rz::real_field::dot_product(const real_field& a)const
{
	return this->inner_product(a);
}

real_field math_rz::real_field::distance(const real_field& a)const
{
	return std::abs(_v-a._v);
}

real_field math_rz::real_field::metric(const real_field& a) const
{
	return std::abs(_v - a._v);
}

real_field math_rz::real_field::inv() const
{
	return 1 / _v;
}

bool real_field::is_zero() const
{
	if(exact) return _v==0;
}

bool real_field::is_one() const
{
	return _v==1;
}


std::ostream& operator<<(std::ostream& H, const real_field& a)
{
	return H << a._v;
}

bool real_field::exact = true;
