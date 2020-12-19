#include "real_field.h"

using namespace math_rz;
//real_field::real_field(double m) :_v(m)
//{
//}

real_field::real_field(long double m):_v(m)
{
}

real_field::real_field(double m):_v(m)
{
}

real_field::real_field(float m):_v(m)
{
}

real_field::real_field(int m):_v(m)
{
}

real_field::real_field(unsigned long long m):_v(m)
{
}

real_field real_field::abs() const
{
	return std::abs(_v);
}

bool real_field::is_zero() const
{
	return _v==0;
}

bool real_field::is_one() const
{
	return _v==1;
}


std::ostream& operator<<(std::ostream& H, const real_field& a)
{
	return H << a._v;
}