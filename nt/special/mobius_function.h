//
// Created by ramizouari on 13/10/2021.
//

#ifndef MATHSLIBRARY_MOBIUS_FUNCTION_H
#define MATHSLIBRARY_MOBIUS_FUNCTION_H
#include "nt/arithmetic_function.h"

namespace math_rz::nt
{
    class mobius_function: public multiplicative_function
    {
    public:
        using multiplicative_function::multiplicative_function;
        integer operator()(const integer &n) const;
    };
    using möbius_function=mobius_function;
}

#endif //MATHSLIBRARY_MOBIUS_FUNCTION_H
