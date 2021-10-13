//
// Created by ramizouari on 07/10/2021.
//

#include "factoriser.h"

namespace math_rz::nt{
    std::vector<integer> factoriser::prime_factors(integer n)
    {
        if(n==1)
            return {};
        auto s = smallest_divisor(n);
        auto P=prime_factors(n/s);
        if(P.empty() || P.back()>s)
            P.push_back(s);
        return P;
    }

    std::unordered_map<integer, integer> factoriser::factorise(integer n)
    {
        std::unordered_map<integer,integer> mapper;
        while(n>1)
        {
            auto s= smallest_divisor(n);
            mapper[s]++;
            n/=s;
        }
        return mapper;
    }

    bool factoriser::is_prime(integer n) {
        return smallest_divisor(n)==n;
    }

    integer factoriser::sum_over_divisors(int n, const std::function<integer(integer)> &f)
    {
        return reduce_over_divisors(n,f,std::plus<integer>());
    }

    integer factoriser::reduce_over_divisors(int n, const std::function<integer(integer)> &f,const std::function<integer(integer,integer)>&O, int o,
                                             const std::vector<integer> &P)
    {
        integer R=f(n);
        int m=P.size();
        for(int i=o;i<m;i++) if(n%P[i]==0)
                R+= reduce_over_divisors(n/P[i],f,O,i,P);
        return R;
    }
    integer factoriser::reduce_over_divisors(int n,const std::function<integer(integer)> &f,const std::function<integer(integer,integer)> &O)
    {
        auto P= prime_factors(n);
        return reduce_over_divisors(n,f,O,0,P);
    }

    integer factoriser::multiply_over_divisors(int n, const std::function<integer(integer)> & f) {
        return reduce_over_divisors(n,f,std::multiplies<integer>());
    }

    integer naive_factoriser::smallest_divisor(integer n)
    {
        auto L=sqrtl(n);
        for(int i=2;i<=L;i++) if(n%i==0)
                return i;
        return n;
    }

    simple_sieve::simple_sieve(integer _N) :N(_N),smallest_d(N+1,0)
    {
        for(int i=2;i<=N;i++)
        {
            integer L=sqrt(i);
            for(auto p:p_list)
            {
                if(i%p==0)
                {
                    smallest_d[i]=p;
                    break;
                }
                if (p > L)
                    break;
            }
            if(smallest_d[i]==0) {
                smallest_d[i] = i;
                p_list.push_back(i);
            }
        }
    }

    integer simple_sieve::smallest_divisor(integer n) {
        return smallest_d[n];
    }

    eratosthenes_sieve::eratosthenes_sieve(integer _N) :N(_N),smallest_d(N+1)
    {
        std::vector<bool> is_prime(N+1,true);
        for(integer i=2;i<=N;i++) if(is_prime[i])
        {
            smallest_d[i]=i;
            for(integer j=i;i*j<=N;j++)
            {
                is_prime[i*j]=false;
                smallest_d[i*j]=i;
            }
        }
    }

    integer eratosthenes_sieve::smallest_divisor(integer n) {
        return smallest_d[n];
    }
}