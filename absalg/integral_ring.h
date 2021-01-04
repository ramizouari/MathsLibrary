#pragma once
#include "ring.h"
#include <utility>
namespace math_rz {
    class integral_ring :
        virtual public ring
    {

    };

    inline std::pair<long long, long long> bezout(const long long& a, const long long& b)
    {
        std::pair<long long, long long> P;
        if (a < b)
        {
            P = bezout(b, a);
            return { P.second,P.first };
        }
        long long r0=a, r1=b,t0=0,t1=1,s0=1,s1=0,w1,w2,w3,q;
        while (r1 != 0)
        {
            w1 = r0;
            w2 = t0;
            w3 = s0;
            r0 = r1;
            s0 = s1;
            t0 = t1;
            q = w1 / r1;
            r1 = w1 -q*r1;
            t1 = w2 - q * t1;
            s1 = w3 - q * s1;
        }
        return { s0,t0 };
    }
}