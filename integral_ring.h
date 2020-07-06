#pragma once
#include "ring.h"
class integral_ring :
    virtual public ring
{
public:
    integral_ring() {}
    virtual integral_ring&& div(const integral_ring&) = 0;
    virtual integral_ring&& mod(const integral_ring&) = 0;
    static integral_ring&& gcd(const integral_ring&, const integral_ring&);
    static integral_ring&& lcm(const integral_ring&, const integral_ring&);
};

