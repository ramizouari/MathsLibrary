#pragma once
#include "ring.h"
#include <utility>
namespace math_rz {
    class integral_ring :
        virtual public ring
    {

    };

    template<typename I>
    std::pair<I, I> bezout(const I& a, const I& b)
    {
        std::pair<I, I> P;
        if (a < b)
        {
            P = bezout(b, a);
            return { P.second,P.first };
        }
        I r0=a, r1=b,t0=0,t1=1,s0=1,s1=0,w1,w2,w3,q;
        while (!r1.is_zero())
        {
            w1 = r0;
            w2 = t0;
            w3 = s0;
            r0 = r1;
            s0 = s1;
            t0 = t1;
            q = w1.div(r1);
            r1 = w1 -q*r1;
            t1 = w2 - q * t1;
            s1 = w3 - q * s1;
        }
        return { s0,t0 };
    }


    template<typename I>
    std::tuple<I,I,I> gcd_bezout(const I& a, const I& b)
    {
        std::tuple<I, I> P;
        if (a < b)
        {
            P = bezout(b, a);
            return { std::get<1>(P),std::get<0>(P),std::get<2>(P) };
        }
        I r0=a, r1=b,t0=0,t1=1,s0=1,s1=0,w1,w2,w3,q;
        while (!r1.is_zero())
        {
            w1 = r0;
            w2 = t0;
            w3 = s0;
            r0 = r1;
            s0 = s1;
            t0 = t1;
            q = w1.div(r1);
            r1 = w1 -q*r1;
            t1 = w2 - q * t1;
            s1 = w3 - q * s1;
        }
        return { s0,t0,s0*a+t0*b };
    }

   
}