//
// Created by ramizouari on 13/10/2021.
//

#include "liouville_function.h"

using namespace math_rz::nt;
integer liouville_function::operator()(const integer &n) const {
    return F.factorise(n).size()%2==0?1:-1;
}
