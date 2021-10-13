//
// Created by ramizouari on 12/10/2021.
//

#ifndef MATHSLIBRARY_ARITHMETIC_FUNCTION_H
#define MATHSLIBRARY_ARITHMETIC_FUNCTION_H
#include "factoriser.h"
#include "analysis/function.h"

namespace math_rz::nt
{
    class arithmetic_function :public math_rz::analysis::abstract_function<integer,integer>
    {
        protected:
            factoriser &F;
    public:
        explicit arithmetic_function(factoriser &_F);
    };

    class multiplicative_function : public arithmetic_function
    {

    public:
        using arithmetic_function::arithmetic_function;
    };
}



#endif //MATHSLIBRARY_ARITHMETIC_FUNCTION_H
