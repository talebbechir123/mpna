// CSRMatrix.hpp
#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#include <vector>

class CSRMatrix {
public:
    int n;                      // dimension
    int nnz;                    // number of non-zeros
    std::vector<int> row_ptr;   // size n+1
    std::vector<int> col_ind;   // size nnz
    std::vector<double> val;    // size nnz

   
    CSRMatrix(int n_, int nnz_)
      : n(n_), nnz(nnz_),
        row_ptr(n_ + 1, 0),
        col_ind(nnz_),
        val(nnz_)
    {}

    
    void multiply(const std::vector<double> &x,
                  std::vector<double> &y) const;
};

#endif
