//
// Created by ramizouari on 13/10/2021.
//

#ifndef MATHSLIBRARY_DIVISOR_FUNCTION_H
#define MATHSLIBRARY_DIVISOR_FUNCTION_H
#include "nt/arithmetic_function.h"

namespace math_rz::nt
{
    class divisor_function:public multiplicative_function
    {
    private:
        integer s;
    public:
        explicit divisor_function(factoriser &_F,integer _s);
        integer operator()(const integer &n) const override;
    };

    class count_divisor_function:public divisor_function
    {
    public:
        explicit count_divisor_function(factoriser &_F);
        integer operator()(const integer &n) const override;
    };
}


#endif //MATHSLIBRARY_DIVISOR_FUNCTION_H
