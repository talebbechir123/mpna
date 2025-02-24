#include "Laplacian2D.h"

// We assume your CSRMatrix has something like:
//    int n;  // dimension
//    std::vector<int> row_ptr;
//    std::vector<int> col_ind;
//    std::vector<double> val;
//    // plus a constructor and multiply() method in CSRMatrix.cpp

CSRMatrix createLaplacian2D(int size)
{
    
    int n = size * size;
    int nnz = 5 * n - 4 * size;

    CSRMatrix A(n, nnz);

    // row_ptr[0] must start at 0
    A.row_ptr[0] = 0;

    int idx = 0;  // track index into col_ind / val

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int row = i*size + j;

            // Center diagonal entry = +4
            A.val[idx]     =  4.0;
            A.col_ind[idx] =  row;
            idx++;

            // Left neighbor
            if (j > 0){
                A.val[idx]     = -1.0;
                A.col_ind[idx] = row - 1;
                idx++;
            }
            // Right neighbor
            if (j < size - 1){
                A.val[idx]     = -1.0;
                A.col_ind[idx] = row + 1;
                idx++;
            }
            // Top neighbor
            if (i > 0){
                A.val[idx]     = -1.0;
                A.col_ind[idx] = row - size;
                idx++;
            }
            // Bottom neighbor
            if (i < size - 1){
                A.val[idx]     = -1.0;
                A.col_ind[idx] = row + size;
                idx++;
            }

            // We have written all non-zeros for this row.
            // The next row's start index is 'idx'
            A.row_ptr[row + 1] = idx;
        }
    }

    // Now A should be correctly filled with an SPD Laplacian matrix.
    return A;
}
