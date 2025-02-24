#include "CSRMatrix.h"
#include <vector>
#include <cmath>
#include <iostream>

/**
 * Solve A*x = b by Conjugate Gradient.
 *
 * @param A         CSRMatrix (assumed SPD)
 * @param b         RHS vector (length A.n)
 * @param x         On input: initial guess; on output: solution
 * @param max_iter  Maximum iterations
 * @param tol       Convergence tolerance (norm of residual)
 */
void conjugateGradientCSR(const CSRMatrix &A,
                          const std::vector<double> &b,
                          std::vector<double> &x,
                          int max_iter,
                          double tol)
{
    int n = A.n;
    if((int)b.size() != n || (int)x.size() != n) {
        std::cerr << "Error: dimension mismatch in CG.\n";
        return;
    }

    // Temporary vectors
    std::vector<double> r(n, 0.0);
    std::vector<double> p(n, 0.0);
    std::vector<double> Ap(n, 0.0);

    // r = b - A*x
    A.multiply(x, r);
    for(int i = 0; i < n; i++){
        r[i] = b[i] - r[i];
    }
    // p = r
    p = r;

    // rs_old = r^T r
    double rs_old = 0.0;
    for(int i = 0; i < n; i++){
        rs_old += r[i]*r[i];
    }

    for(int k = 0; k < max_iter; k++){
        // Ap = A*p
        A.multiply(p, Ap);

        // alpha = (r^T r) / (p^T Ap)
        double pAp = 0.0;
        for(int i = 0; i < n; i++){
            pAp += p[i]*Ap[i];
        }
        double alpha = rs_old / pAp;

        // x = x + alpha*p
        for(int i = 0; i < n; i++){
            x[i] += alpha * p[i];
        }

        // r = r - alpha*Ap
        for(int i = 0; i < n; i++){
            r[i] -= alpha * Ap[i];
        }

        // rs_new = r^T r
        double rs_new = 0.0;
        for(int i = 0; i < n; i++){
            rs_new += r[i]*r[i];
        }

        // Check convergence
        if(std::sqrt(rs_new) < tol){
            //std::cout << "CG converged in " << k+1 << " iterations.\n";
            break;
        }

        // beta = rs_new / rs_old
        double beta = rs_new / rs_old;
        // p = r + beta*p
        for(int i = 0; i < n; i++){
            p[i] = r[i] + beta * p[i];
        }

        rs_old = rs_new;
    }
}
