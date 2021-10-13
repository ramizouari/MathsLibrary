//
// Created by ramizouari on 13/10/2021.
//

#include "divisor_function.h"

using namespace math_rz::nt;

math_rz::integer divisor_function::operator()(const integer &n) const
{
    integer R=1;
    auto mapper=F.factorise(n);
    if(s==0) for(auto [_,m]:mapper)
        R*=(m+1);
    else for(auto [p,m]:mapper)
        R*=(pow(p,(m+1)*s)-1)/(pow(p,s)-1);
    return R;
}

divisor_function::divisor_function(factoriser &_F, math_rz::integer _s): multiplicative_function(_F),s(_s) {

}

count_divisor_function::count_divisor_function(factoriser &_F) : divisor_function(_F, 0) {

}

math_rz::integer count_divisor_function::operator()(const math_rz::integer &n) const {
    integer R=1;
    for(auto [_,m]:F.factorise(n))
        R*=(m+1);
    return R;
}
