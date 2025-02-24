#ifndef LAPLACIAN_2D_H
#define LAPLACIAN_2D_H

#include "CSRMatrix.h"

/**
 * Build a 2D Laplacian matrix (5-point stencil) of dimension (size*size),
 * stored in CSR format.
 *
 * @param size  The number of points in one dimension of a 2D grid.
 *              The total system size will be size^2.
 * @return      A CSRMatrix representing the 5-pt Laplacian.
 */
CSRMatrix createLaplacian2D(int size);

#endif // LAPLACIAN_2D_HPP
