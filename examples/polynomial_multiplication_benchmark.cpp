//
// Created by ramizouari on 08/07/2021.
//

#include "linalg/matrix.h"
#include "poly/multiplicator/multiplicator.h"
#include "poly/polynomial.h"
#include "poly/structure/inner_product.h"
#include "complex.h"
#include "absalg/cyclic.h"


using namespace math_rz;
using namespace math_rz::linalg;
using namespace math_rz::analysis;
using K = math_rz::complex;
using M = math_rz::linalg::matrix<K, 3,5>;
constexpr int dimension = 3;
constexpr int points=5e4;

int main()
{
    using namespace std::complex_literals;
    std::vector<real_field> S(points);
    for (int i = 0; i < points; i++)
        S[i] = 1;
    poly::polynomial<real_field> p(S);
    std::chrono::time_point t1=std::chrono::system_clock::now();
    auto R1= p * p;
    decltype(p)::set_multiplicator(new poly::multiplicator::fast_multiplicator<real_field>());
    decltype(p)::set_structure(new poly::structure::L2_vect_inner_product<real_field>());
    std::chrono::time_point t2=std::chrono::system_clock::now();
    std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count())/1000 << std::endl;
    auto R2 = p * p;
    std::chrono::time_point t3=std::chrono::system_clock::now();
    std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count())/1000 << std::endl;
    std::cout << (R2 - R1).norm() << "\n" << R1;
    return false;
}