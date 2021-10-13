//
// Created by ramizouari on 12/10/2021.
//

#include "totient.h"
#include "absalg/ring.h"
using namespace math_rz::nt;

math_rz::integer euler_totient::operator()(const integer &a) const
{
    auto mapper=F.factorise(a);
    integer R=1;
    for(auto [p,m]:mapper)
        R*=pow(p,m-1)*(p-1);
    return R;
}

math_rz::integer carmichael_totient::operator()(const math_rz::integer &a) const
{
    auto mapper=F.factorise(a);
    integer R=1;
    for(auto [p,m]:mapper)
        if(p==2 && m>2)
            R=std::lcm<long long>(R,pow(p,m-2)*(p-1));
        else R=std::lcm<long long>(R,pow(p,m-1)*(p-1));
    return R;
}
