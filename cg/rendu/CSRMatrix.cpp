// CSRMatrix.cpp
#include "CSRMatrix.h"
#include <iostream>

void CSRMatrix::multiply(const std::vector<double> &x,
                         std::vector<double> &y) const
{
    if((int)x.size() != n || (int)y.size() != n){
        std::cerr << "Error: Dimension mismatch in CSRMatrix::multiply.\n";
        return;
    }

    for(int i = 0; i < n; i++){
        double sum = 0.0;
        for(int idx = row_ptr[i]; idx < row_ptr[i+1]; idx++){
            sum += val[idx] * x[col_ind[idx]];
        }
        y[i] = sum;
    }
}
