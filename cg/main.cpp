#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

// Declarations for CG, Jacobi, Gauss-Seidel
void conjugateGradientCSR(const class CSRMatrix &A,
                          const std::vector<double> &b,
                          std::vector<double> &x,
                          int max_iter,
                          double tol);

#include "iterative_methods.h"
#include "CSRMatrix.h"
#include "CSRMatrix.h"
#include "iterative_methods.h"
#include "Laplacian2D.h"  


int main() {
    // Suppose we want a Laplacian for a 10x10 grid => dimension = 100
    int size = 10;

    // Build the matrix
    CSRMatrix A = createLaplacian2D(size);

    // Build the RHS vector b, for example, all ones:
    std::vector<double> b(A.n, 1.0);

    // Prepare an initial guess x = zero
    std::vector<double> x(A.n, 0.0);

    // 1) Test Conjugate Gradient
    int max_iter = 1000;
    double tol = 1e-10;

    auto start = std::chrono::high_resolution_clock::now();
    conjugateGradientCSR(A, b, x, max_iter, tol /*, etc.*/);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur_cg = end - start;

    std::cout << "CG finished in " << dur_cg.count() << " seconds.\n";

    // 2) Reset x, test Jacobi
    std::fill(x.begin(), x.end(), 0.0);

    start = std::chrono::high_resolution_clock::now();
    jacobiCSR(A, b, x, max_iter, tol /*, etc.*/);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur_jacobi = end - start;

    std::cout << "Jacobi finished in " << dur_jacobi.count() << " seconds.\n";

    // 3) Reset x, test Gauss-Seidel
    std::fill(x.begin(), x.end(), 0.0);

    start = std::chrono::high_resolution_clock::now();
    gaussSeidelCSR(A, b, x, max_iter, tol /*, etc.*/);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur_gs = end - start;

    std::cout << "Gauss-Seidel finished in " << dur_gs.count() << " seconds.\n";

    // Possibly print results or residual norms, etc.
    return 0;
}
