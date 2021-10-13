//
// Created by ramizouari on 17/09/2021.
//

#ifndef MATHSLIBRARY_BIG_INT_H
#define MATHSLIBRARY_BIG_INT_H
#include <array>
namespace math_rz{
    template<std::integral I, int d>
    class big_int {
        inline static constexpr int n=8*sizeof(I),p=(d+n-1)/n;
        std::array<I,p> v;
    public:

        big_int(I a=0)
        {
            v[0]=a;
        }
        template<int s>
        big_int(big_int<I,s>A)
        {
            int m=std::min(d,s);
            for(int i=0;i<m;i++)
                v[i]=A.v[i];
        }
        template<int s>
        big_int operator+=(big_int<I,s> A)
        {
            int m=std::min(d,s);
            I C=0;
            for(int i=0;i<m;i++)
            {
                if(C)
                    v[i]+=C;
                C=C&&v[i]==0;
                v[i]+=A.v[i];
                C=C||(v[i]<A.v[i]);
            }
            return *this;
        }

        big_int operator+=(I A)
        {
            return *this+=big_int<I,1>(A);
        }
    };
}


#endif //MATHSLIBRARY_BIG_INT_H
