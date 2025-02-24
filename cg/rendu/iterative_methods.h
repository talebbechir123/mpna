#ifndef ITERATIVE_METHODS_H
#define ITERATIVE_METHODS_H

#include "CSRMatrix.h"
#include <vector>

/**
 * Solve A*x = b by Jacobi iteration (CSR format).
 *
 * @param A         CSRMatrix (assumed SPD or at least diagonally dominant)
 * @param b         RHS vector
 * @param x         On input: initial guess; on output: solution
 * @param max_iter  Maximum iterations
 * @param tol       Convergence tolerance
 */
void jacobiCSR(const CSRMatrix &A,
               const std::vector<double> &b,
               std::vector<double> &x,
               int max_iter,
               double tol);

/**
 * Solve A*x = b by Gauss-Seidel iteration (CSR format).
 *
 * @param A         CSRMatrix
 * @param b         RHS vector
 * @param x         On input: initial guess; on output: solution
 * @param max_iter  Maximum iterations
 * @param tol       Convergence tolerance
 */
void gaussSeidelCSR(const CSRMatrix &A,
                    const std::vector<double> &b,
                    std::vector<double> &x,
                    int max_iter,
                    double tol);

#endif
