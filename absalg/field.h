#pragma once
#include "integral_ring.h"
#include <concepts>

namespace math_rz
{
    class field :
        public integral_ring
    {
    };
    namespace field_constraints
    {

        template<typename F>
        concept field = requires(const F & a, const F&b)
        {
            a + b;
            a - b;
            a* b;
            a / b;
            -a;
            a.inv();
        };

        template<typename F>
        concept has_abs = field<F> && requires(F a)
        {
            a.abs();
        };
        
        /*
        * An absolute value of a field K is a function abs:K->reals such that
        * - abs is positive definite
        * - abs is multiplicative
        * - abs verifies the triangle-inequality
        * Such field is fundamental for the definition of the notion of a norm
        * If defined, it is expected that a norm is compatible with the absolute value
        */
        template<typename K>
        concept field_with_abs = field<K> && requires(const K & a)
        {
            a.abs();
            a.norm();
        };
        
        /*
        * We define a field with conjugate as a field with an involutive field automorphism L:
        * - The real numbers with L=id
        * - The complex numbers with L the natural conjugate
        * If the field has an absolute value, it is expected that the conjugate is compatible with the absolute value:
        * - |w|² = w * w.conj()
        */
        template<typename F>
        concept field_with_conj = field<F>&&requires(F a)
        {
            a.conj();
        };

    }
}