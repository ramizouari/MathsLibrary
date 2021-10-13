//
// Created by ramizouari on 08/10/2021.
//


#include "utils.h"
#include "absalg/integral_ring.h"
#include <stack>

namespace math_rz::nt
{
    integer chinese_remainder(const std::vector<std::pair<integer,integer>> &S)
    {
        std::stack<std::pair<integer,integer>> Q;
        for(auto s:S)
            Q.push(s);
        while(Q.size() > 1)
        {
            auto [a1,p1]=Q.top();
            Q.pop();
            auto [a2,p2]=Q.top();
            Q.pop();
            auto [k1,k2]=bezout(p1,p2);
            k2*=(a1-a2);
            Q.push({(k2*p2+a2)%(p1*p2),p1*p2});
        }
        return Q.top().first;
    }
    integer chinese_remainder(std::vector<integer> A,std::vector<integer> P)
    {
        std::vector<std::pair<integer,integer>> S;
        int n=A.size(),m=P.size();
        S.reserve(n);
        for(int i=0;i<n;i++)
            S.emplace_back(A[i],P[i]);
        return chinese_remainder(S);
    }

}