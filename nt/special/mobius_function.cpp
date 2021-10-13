//
// Created by ramizouari on 13/10/2021.
//

#include "mobius_function.h"

math_rz::integer math_rz::nt::mobius_function::operator()(const math_rz::integer &n) const
{
    auto mapper=F.factorise(n);
    for(auto [_,m]:mapper)
        if(m>1)
            return 0;
    return mapper.size()%2==0?1:-1;
}
