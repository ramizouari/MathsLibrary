//
// Created by ramizouari on 07/10/2021.
//

#ifndef MATHSLIBRARY_FACTORISER_H
#define MATHSLIBRARY_FACTORISER_H
#include "integer.h"
#include <cmath>
#include <unordered_map>
#include <functional>

namespace math_rz::nt
{
    class factoriser
    {
        integer reduce_over_divisors(int n,const std::function<integer(integer)> &f,const std::function<integer(integer,integer)>&O,
                                     int o,const std::vector<integer> &P);
    public:
        virtual integer smallest_divisor(integer n) = 0;
        virtual bool is_prime(integer n);
        virtual std::vector<integer> prime_factors(integer n);
        virtual std::unordered_map<integer,integer> factorise(integer n);
        integer sum_over_divisors(int n,const std::function<integer(integer)> &f);
        integer multiply_over_divisors(int n,const std::function<integer(integer)> &f);
        integer reduce_over_divisors(int n,const std::function<integer(integer)> &f,const std::function<integer(integer,integer)> &O);
    };

    class naive_factoriser:public factoriser
    {
        integer smallest_divisor(integer n) override;
    };

    class simple_sieve:public factoriser
    {
        integer N;
        std::vector<integer> p_list,smallest_d;
    public:
        explicit simple_sieve(integer _N);
        integer smallest_divisor(integer n) override;
    };

    class eratosthenes_sieve:public factoriser
    {
        integer N;
        std::vector<integer> p_list,smallest_d;
    public:
        eratosthenes_sieve(integer _N);
        integer smallest_divisor(integer n) override;

    };
}


#endif //MATHSLIBRARY_FACTORISER_H
