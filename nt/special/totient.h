//
// Created by ramizouari on 12/10/2021.
//

#ifndef MATHSLIBRARY_TOTIENT_H
#define MATHSLIBRARY_TOTIENT_H

#include "nt/factoriser.h"
#include "nt/arithmetic_function.h"
namespace math_rz::nt
{
    class euler_totient:public multiplicative_function
    {
    public:
        using multiplicative_function::multiplicative_function;
        integer operator()(const integer& a) const override;
    };

    class carmichael_totient:public multiplicative_function
    {
    public:
        using multiplicative_function::multiplicative_function;
        integer operator()(const integer& a) const override;
    };
}


#endif //MATHSLIBRARY_TOTIENT_H
