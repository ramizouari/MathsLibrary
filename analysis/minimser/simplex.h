#pragma once
#include "linalg/finite_dimensional_vector_space.h"
#include "linalg/matrix.h"
namespace math_rz::analysis
{
    template<int n, int m>
    linalg::finite_dimensional_vector_space<real_field,n> arg_max(
        const linalg::finite_dimensional_vector_space<real_field,n>& _Z, 
        const linalg::matrix<real_field,n, m>& _A, 
        linalg::finite_dimensional_vector_space<real_field,m> b)
    {
        linalg::finite_dimensional_vector_space < real_field, n + m > Z;
        for (int i = 0; i < n; i++)
            Z[i] = _Z[i];
        linalg::matrix<real_field,n, n + m>A;
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
            A[i][j] = _A[i][j];
        for (int i = m; i < n + m; i++)
            A[i - m][i] = 1;
        linalg::finite_dimensional_vector_space<real_feidl,n>U;
        while (Z.max() > 0)
        {
            auto q = Z.arg_max();
            int p = -1;
            auto c = b.point_wise_divide(A.get_column(q));
            for (int k = 0; k < n; k++)
                if (A[k][q] > 0 && c[k] >= 0 && (p == -1 || c[k] < c[p]))
                    p = k;
            if (p == -1)
                break;
            for (int i = 0; i < n; i++) if (i != p)
            {
                b[i] -= (A[i][q] / A[p][q]) * b[p];
                A[i] -= (A[i][q] / A[p][q]) * A[p];
                A[i][q] = 0;

            }
            Z -= (Z[q] / A[p][q]) * A[p];
            Z[q] = 0;
        }
        linalg::matrix<real_field,n, n> P;
        linalg::finite_dimensional_vector_space<real_field,n + m> h;
        for (int i = 0; i < m; i++)
            h[i] = b[i];
        linalg::matrix<real_field,n + m, n + m> Q;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n + m; j++)
                Q[i][j] = A[i][j];
        int r = 0;
        for (int i = 0; i < n + m ; i++)
            if (Z[i] < 0)
            {
                Q[n + r++][i] = 1;
                for (int j = 0; j < n; j++)
                    A[j][i] = 0;
            }
        linalg::finite_dimensional_vector_space<real_field,n> d;
        auto g = Q.solve(h);
        for (int i = 0; i < n; i++)
            d[i] = g[i];
        return d;
    }
}