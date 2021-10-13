//
// Created by ramizouari on 13/10/2021.
//

#ifndef MATHSLIBRARY_LIOUVILLE_FUNCTION_H
#define MATHSLIBRARY_LIOUVILLE_FUNCTION_H
#include "nt/arithmetic_function.h"
namespace math_rz::nt
{
    class liouville_function :public multiplicative_function
    {
        public:
        usign multiplicative_function::multiplicative_function;
        integer operator()(const integer &n) const override;
    };
}


#endif //MATHSLIBRARY_LIOUVILLE_FUNCTION_H
