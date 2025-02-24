
#ifndef GMRES_H
#define GMRES_H


#include <vector>
#include <cmath>
#include <iostream>
#include "CSRMatrix.h"
#include "Conjugate-gradient.h"

// making a GMRES solver

template<typename IdType, typename ScalarType>
int gmres(int n, const CSRMatrix<IdType, ScalarType>& A, const std::vector<ScalarType>& b, std::vector<ScalarType>& x, ScalarType tol, int max_iter, int m)
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

        // check for convergence
        ScalarType r_norm   = 0.0;  // residual norm    
        ScalarType b_norm   = 0.0;  // rhs norm
        ScalarType rel_norm = 0.0;  // relative norm
        for (int i = 0; i < n; ++i)
        {
            r_norm += r[i] * r[i];
            b_norm += b[i] * b[i];
        }
        r_norm = std::sqrt(r_norm);
        b_norm = std::sqrt(b_norm);
        rel_norm = r_norm / b_norm;
        if (rel_norm < tol)
        {
            return iter;
        }

        // solve the least squares problem
        std::vector<std::vector<ScalarType>> H(m + 1, std::vector<ScalarType>(m, 0.0));
        std::vector<ScalarType> c(m + 1, 0.0);
        std::vector<ScalarType> s(m + 1, 0.0);
        std::vector<ScalarType> y(m + 1, 0.0);
        std::vector<ScalarType> e1(m + 1, 0.0);
        e1[0] = 1.0;

        for (int i = 0; i < m; ++i)
        {
            // construct the i-th column of the Hessenberg matrix
            for (int j = 0; j < n; ++j)
            {
                H[i][j] = 0.0;
            }
            for (int j = 0; j < n; ++j)
            {
                H[i][j] = p[j] * r[j];
            }

            // apply Givens rotation
            for (int k = 0; k < i; ++k)
            {
                ScalarType temp = c[k] * H[k][i] + s[k] * H[k + 1][i];
                H[k + 1][i] = -s[k] * H[k][i] + c[k] * H[k + 1][i];
                H[k][i] = temp;
            }

            // construct the i-th Givens rotation
            ScalarType r1 = std::sqrt(H[i][i] * H[i][i] + H[i + 1][i] * H[i + 1][i]);
            c[i] = H[i][i] / r1;
            s[i] = -H[i + 1][i] / r1;

            // apply the i-th Givens rotation
            H[i][i] = c[i] * H[i][i] + s[i] * H[i + 1][i];
            H[i + 1][i] = 0.0;

            // apply the i-th Givens rotation to the right-hand side
            ScalarType temp = c[i] * r[i] + s[i] * r[i + 1];
            r[i + 1] = -s[i] * r[i] + c[i] * r[i + 1];
            r[i] = temp;

            // apply the i-th Givens rotation to the right-hand side
            temp = c[i] * b[i] + s[i] * b[i + 1];
            b[i + 1] = -s[i] * b[i] + c[i] * b[i + 1];
            b[i] = temp;
        }

        // solve the least squares problem
        for (int i = m - 1; i >= 0; --i)
        {
            y[i] = b[i];
            for (int j = i + 1; j < m; ++j)
            {
                y[i] -= H[i][j] * y[j];
            }
            y[i] /= H[i][i];
        }

        // update the solution
        for (int i = 0; i < n; ++i)
        {
            x[i] += y[i] * p[i];
        }

        // update the residual
        for (int i = 0; i < n; ++i)
        {
            r[i] -= y[i] * Ap[i];
        }

        // check for convergence
        r_norm = 0.0;
        for (int i = 0; i < n; ++i)
        {
            r_norm += r[i] * r[i];
        }

        r_norm = std::sqrt(r_norm);
        rel_norm = r_norm / b_norm;
        if (rel_norm < tol)
        {
            return iter;
        }
    }

    return max_iter;

}

#endif // GMRES_H
