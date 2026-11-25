# Advanced Numerical Programming — Dense, Sparse and Iterative Linear Algebra

My solutions to the lab series of the *Méthodes de Programmation Numérique Avancée* (MPNA) course, Master CHPS: dense matrix products, sparse CSR kernels, classical and Krylov iterative solvers, and an eigenvalue application. Written in C and C++.

Lab statements and course framework by Cédric Chevalier (BSD 2-Clause, see [`LICENSE`](LICENSE)); lecture slides are in [`slides/`](slides/).

## Contents

| Folder | Topic | What is implemented |
|---|---|---|
| [`01-matrix-matrix/`](01-matrix-matrix/) | Dense matrix multiplication (C) | Naive and loop-reordered implementations with timing over matrix sizes |
| [`02-iterative-relaxation/`](02-iterative-relaxation/) | Iterative solvers (C++20, header-only) | Templated CSR matrix, 2D Laplacian generator, **Jacobi**, **Gauss-Seidel**, **Conjugate Gradient**, restarted **GMRES(m)** |
| [`03-sparse-matvec/`](03-sparse-matvec/) | Sparse matrix-vector product (C) | CSR SpMV reading Matrix Market files, **power iteration** for the dominant eigenvalue; MPI version started in `mat-vec-dist.c` |
| [`cg/`](cg/) | Conjugate Gradient on a 2D Poisson problem (C++17) | CSR solver for the 2D Laplacian, comparison with Jacobi and Gauss-Seidel, analysis notebook `CG2D.ipynb` |

## Build

Labs 01 and 02 use CMake from the repository root:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/01-matrix-matrix/matrix-matrix <n> <m> <p>
./build/02-iterative-relaxation/iter_laplacian
```

The other labs have their own Makefile:

```bash
cd 03-sparse-matvec && make && ./mat-vec   # reads data/bcsstk03/bcsstk03.mtx (path set in main)
cd cg && make && ./solver
```

Test matrices in Matrix Market format are available from the [SuiteSparse Matrix Collection](https://sparse.tamu.edu/).

## Concepts covered

- Cache-aware loop ordering and data layout for dense kernels
- Compressed Sparse Row (CSR) storage and SpMV
- Convergence behaviour of stationary methods (Jacobi, Gauss-Seidel) versus Krylov methods (CG, GMRES)
- Eigenvalue estimation by power iteration

## Author

**Ahmed Taleb Bechir**
[atbechir.com](https://atbechir.com) · [LinkedIn](https://www.linkedin.com/in/ahmed-taleb-bechir-97a7301a4) · [GitHub](https://github.com/talebbechir123)
