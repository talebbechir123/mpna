#include "iterative_methods.h"
#include <cmath>
#include <vector>
#include <iostream>

void jacobiCSR(const CSRMatrix &A,
               const std::vector<double> &b,
               std::vector<double> &x,
               int max_iter,
               double tol)
{
    int n = A.n;
    if((int)b.size() != n || (int)x.size() != n) {
        std::cerr << "Error: dimension mismatch in jacobiCSR.\n";
        return;
    }

    std::vector<double> x_old(n, 0.0);

    for(int iter = 0; iter < max_iter; iter++){
        x_old = x;  // Keep a copy of the old iterate

        // Update each row i
        for(int i = 0; i < n; i++){
            double diag = 0.0;
            double sum  = 0.0;
            // Row i: loop over non-zeros
            for(int idx = A.row_ptr[i]; idx < A.row_ptr[i+1]; idx++){
                int col = A.col_ind[idx];
                if(col == i) {
                    diag = A.val[idx];
                } else {
                    sum += A.val[idx] * x_old[col];
                }
            }
            x[i] = (b[i] - sum) / diag;
        }

        // Check convergence via difference norm
        double diff_norm2 = 0.0;
        for(int i = 0; i < n; i++){
            double d = x[i] - x_old[i];
            diff_norm2 += d*d;
        }
        if(std::sqrt(diff_norm2) < tol) {
            break;
        }
    }
}

void gaussSeidelCSR(const CSRMatrix &A,
                    const std::vector<double> &b,
                    std::vector<double> &x,
                    int max_iter,
                    double tol)
{
    int n = A.n;
    if((int)b.size() != n || (int)x.size() != n) {
        std::cerr << "Error: dimension mismatch in gaussSeidelCSR.\n";
        return;
    }

    std::vector<double> x_old(n, 0.0);

    for(int iter = 0; iter < max_iter; iter++){
        x_old = x;  // for referencing old values in the same iteration

        // In-place update: row by row
        for(int i = 0; i < n; i++){
            double diag = 0.0;
            double sum  = 0.0;
            for(int idx = A.row_ptr[i]; idx < A.row_ptr[i+1]; idx++){
                int col = A.col_ind[idx];
                if(col == i) {
                    diag = A.val[idx];
                } else {
                    // For col < i, use the newly updated x
                    // For col > i, use the old x
                    sum += A.val[idx] * ((col < i) ? x[col] : x_old[col]);
                }
            }
            x[i] = (b[i] - sum) / diag;
        }

        // Check convergence via difference norm
        double diff_norm2 = 0.0;
        for(int i = 0; i < n; i++){
            double d = x[i] - x_old[i];
            diff_norm2 += d*d;
        }
        if(std::sqrt(diff_norm2) < tol) {
            break;
        }
    }
}
