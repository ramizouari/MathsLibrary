#pragma once
#include "integral_ring.h"
#include <concepts>
class field :
    public integral_ring
{
};
namespace field_constraints
{

    template<typename F>
    concept is_field = std::is_base_of_v<field, F>;
    template<typename F>
    concept has_abs = is_field<F> && requires(F a)
    {
        a.abs();
    } ;

    template<typename F>
    concept is_complex = requires(F a)
    {
        a.conj();
    };

}