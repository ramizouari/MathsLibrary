//
// Created by ramizouari on 12/10/2021.
//

#include "nt/factoriser.h"
#include "nt/special/totient.h"
#include "nt/special/divisor_function.h"
#include "nt/special/mobius_function.h"
#include <iostream>


int main()
{
    using namespace math_rz;
    integer N,n;
    std::cin >> N >> n;
    nt::simple_sieve F(N);
    nt::euler_totient T(F);
    nt::carmichael_totient C(F);
    nt::divisor_function D(F,2);
    nt::möbius_function µ(F);
    std::cout << F.sum_over_divisors(50,[&µ,&D](auto s){return µ(s)*D(50/s);});
}