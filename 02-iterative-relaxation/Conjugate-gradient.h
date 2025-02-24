
#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include <vector>
#include <cmath>
#include <iostream>
#include "CSRMatrix.h"
// making a conjugate gradient solver

template<typename IdType, typename ScalarType>
int conjugate_gradient(int n, const CSRMatrix<IdType, ScalarType>& A, const std::vector<ScalarType>& b, std::vector<ScalarType>& x, ScalarType tol, int max_iter)
{
    std::vector<ScalarType> r(n, 0.0);
    std::vector<ScalarType> p(n, 0.0);
    std::vector<ScalarType> Ap(n, 0.0);

    // r = b - A*x
    for (int i = 0; i < n; ++i)
    {
        ScalarType sum = 0.0;
        for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j)
        {
            sum += A.values[j] * x[A.col_idx[j]];
        }
        r[i] = b[i] - sum;
        p[i] = r[i];
    }

    for (int iter = 0; iter < max_iter; ++iter)
    {
        // Ap = A*p
        for (int i = 0; i < n; ++i)
        {
            ScalarType sum = 0.0;
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j)
            {
                sum += A.values[j] * p[A.col_idx[j]];
            }
            Ap[i] = sum;
        }

        // alpha = (r^T * r) / (p^T * Ap)
        ScalarType alpha = 0.0;
        ScalarType r_dot_r = 0.0;
        ScalarType p_dot_Ap = 0.0;
        for (int i = 0; i < n; ++i)
        {
            r_dot_r += r[i] * r[i];
            p_dot_Ap += p[i] * Ap[i];
        }
        alpha = r_dot_r / p_dot_Ap;

        // x = x + alpha * p
        for (int i = 0; i < n; ++i)
        {
            x[i] += alpha * p[i];
        }

        // r = r - alpha * Ap
        for (int i = 0; i < n; ++i)
        {
            r[i] -= alpha * Ap[i];
        }

        // beta = (r^T * r) / ( r_old^T * r_old)
        ScalarType beta = 0.0;
        ScalarType r_dot_r_new = 0.0;
        for (int i = 0; i < n; ++i)
        {
            r_dot_r_new += r[i] * r[i];
        }
        beta = r_dot_r_new / r_dot_r;
        

        // p = r + beta * p
        for (int i = 0; i < n; ++i)
        {
            p[i] = r[i] + beta * p[i];
        }

        // Simplest stopping criterion: check the norm of the residual
        ScalarType norm = 0.0;
        for (int i = 0; i < n; ++i)
        {
            norm += r[i] * r[i];
        }
        norm = std::sqrt(norm);

        if (norm < tol)
        {
            return iter;
        }
    }
    return max_iter;
}



#endif // CONJUGATE_GRADIENT_H