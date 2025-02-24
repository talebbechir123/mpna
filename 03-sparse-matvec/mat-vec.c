#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mmio.h"
#include <math.h>





//implemeting a csr matrix
typedef struct {
    int n;
    int m;
    int *row_ptr;
    int *col_idx;
    double *val;
} csr_matrix;

// //function to read Matrix Market format file
// void readMMF(char *filename, csr_matrix *A){

//     MM_typecode matcode;
//     FILE *f = fopen(filename,"r");
//     if(f == NULL){
//         printf("Error opening file %s\n",filename);
//         exit(1);
//     }
//     int M, N, nz;   
//     int i, *I, *J;
//     double *val;

//     if(mm_read_banner(f, &matcode) != 0){
//         printf("Could not process Matrix Market banner.\n");
//         exit(1);
//     }

//     if(mm_is_complex(matcode)){
//         printf("Sorry, this application does not support complex numbers.\n");
//         exit(1);
//     }

//     if(mm_is_sparse(matcode)){
//         if(mm_read_mtx_crd_size(f, &M, &N, &nz) != 0){
//             printf("Error reading matrix size.\n");
//             exit(1);
//         }
//     }else{
//         printf("Not a sparse matrix.\n");
//         exit(1);
//     }

//     I = (int *)malloc(nz * sizeof(int));
//     J = (int *)malloc(nz * sizeof(int));
//     val = (double *)malloc(nz * sizeof(double));
    
//     for(i = 0; i < nz; i++){
//         fscanf(f,"%d %d %lf\n",&I[i],&J[i],&val[i]);
//         I[i]--;  //adjust from 1-based to 0-based
//         J[i]--;
//     }

//     fclose(f);

//     A->n = M;
//     A->m = N;
//     A->row_ptr = (int *)malloc((M+1) * sizeof(int));
//     A->col_idx = (int *)malloc(nz * sizeof(int));
//     A->val = (double *)malloc(nz * sizeof(double));

//     for(i = 0; i < M; i++){
//         A->row_ptr[i] = 0;
//     }

//     for(i = 0; i < nz; i++){
//         A->row_ptr[I[i]]++;
//     }

//     for(i = 0; i < M; i++){
//         A->row_ptr[i+1] += A->row_ptr[i];
//     }

//     for(i = 0; i < nz; i++){
//         int row = I[i];
//         int idx = A->row_ptr[row];
//         A->col_idx[idx] = J[i];
//         A->val[idx] = val[i];
//         A->row_ptr[row]++;
//     }

//     int last = 0;
//     for(i = 0; i <= M; i++){
//         int temp = A->row_ptr[i];
//         A->row_ptr[i] = last;
//         last = temp;
//     }

//     //matrix is now in csr format
//     printf("Matrix read successfully\n");

//     free(I);
//     free(J);
//     free(val);

    
// }



void readMMF(const char *filename, csr_matrix *A) {

    MM_typecode matcode;
    FILE *f = fopen(filename, "r");
    if (!f) { /* handle error */ }

    // Read banner
    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    // Check if it's sparse + real + symmetric
    if (!mm_is_sparse(matcode) || !mm_is_real(matcode) || !mm_is_symmetric(matcode)) {
        printf("Matrix must be sparse, real, and symmetric.\n");
        exit(1);
    }

    // Read dimensions and number of entries in the triangular part
    int M, N, nz;
    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
        printf("Error reading matrix size.\n");
        exit(1);
    }

    // Allocate temporary COO arrays (for the triangular portion)
    int *I = (int*)malloc(nz * sizeof(int));
    int *J = (int*)malloc(nz * sizeof(int));
    double *val = (double*)malloc(nz * sizeof(double));

    // Read the (i, j, value) lines
    for (int i = 0; i < nz; i++) {
        fscanf(f, "%d %d %lf", &I[i], &J[i], &val[i]);
        // convert to 0-based
        I[i]--;
        J[i]--;
    }
    fclose(f);

    // Count how many off-diagonal entries
    int offDiag = 0;
    for (int k = 0; k < nz; k++) {
        if (I[k] != J[k]) {
            offDiag++;
        }
    }
    int totalEntries = nz + offDiag; // each off-diagonal is mirrored

    // Allocate CSR structure
    A->n = M;
    A->m = N;
    A->row_ptr = (int*)calloc(M+1, sizeof(int));
    A->col_idx = (int*)malloc(totalEntries * sizeof(int));
    A->val     = (double*)malloc(totalEntries * sizeof(double));

    // 1) Count row sizes
    for (int k = 0; k < nz; k++) {
        A->row_ptr[I[k]]++;
        if (I[k] != J[k]) {
            A->row_ptr[J[k]]++;
        }
    }

    // 2) Prefix sum
    for (int row = 0; row < M; row++) {
        A->row_ptr[row+1] += A->row_ptr[row];
    }

    // 3) Fill col_idx and val
    int *row_pos = (int*)calloc(M, sizeof(int));

    for (int k = 0; k < nz; k++) {
        int r = I[k];
        int c = J[k];
        double v = val[k];

        // Put (r, c)
        int destPos = A->row_ptr[r] + row_pos[r];
        A->col_idx[destPos] = c;
        A->val[destPos]     = v;
        row_pos[r]++;

        // If off-diagonal, also mirror
        if (r != c) {
            destPos = A->row_ptr[c] + row_pos[c];
            A->col_idx[destPos] = r;
            A->val[destPos]     = v;
            row_pos[c]++;
        }
    }

    free(row_pos);
    free(I);
    free(J);
    free(val);

    printf("Matrix read successfully in full symmetric form.\n");
}


//function to read a vector from a file
void readVec(char *filename ,double vec[], int n){
    FILE *f = fopen(filename,"r");
    if(f == NULL){
        printf("Error opening file %s\n",filename);
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fscanf(f,"%lf",&vec[i]);
    }
    fclose(f);
}


// Compute the Euclidean norm (L2 norm) of a vector of length n.
double vector_norm(int n, double* x) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * x[i];
    }
    return sqrt(sum);
}

// Scale a vector x (of length n) by a constant factor.
void scale_vector(int n, double* x, double scale) {
    for (int i = 0; i < n; i++) {
        x[i] *= scale;
    }
}

// Compute the dot product of two vectors of length n.
double dot_product(int n, double* x, double* y) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}



//function to multiply a csr matrix with a vector
void matVec(csr_matrix *A, double x[], double y[]){
    int n = A->n;
    for(int i = 0; i < n; i++){
        y[i] = 0;
        for(int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++){
            y[i] += A->val[j] * x[A->col_idx[j]];
        }
    }
}

//-------------------------------------
// Power Iteration to Estimate Largest Eigenvalue
//-------------------------------------
// We use the norm of A*x (with x normalized) to approximate λ.
double PowerIteration(csr_matrix *A, double x[], int max_iter) {
    int n = A->n;
    double *y = (double *)malloc(n * sizeof(double));
    if (y == NULL) {
        fprintf(stderr, "Memory allocation error.\n");
        exit(EXIT_FAILURE);
    }
    
    double lambda = 0.0;
    
    // Initialize x randomly and normalize it.
    srand((unsigned)time(NULL));
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        x[i] = (double)rand() / (double)RAND_MAX;
       
    }
   
    
    // Power iteration loop.
    for (int iter = 0; iter < max_iter; iter++){
        // y = A*x
        matVec(A, x, y);
        
        // Compute norm of y.
        double y_norm = vector_norm(n, y);
        
        // Approximate eigenvalue: since x is normalized,
        // we have A*x ≈ λ*x so that ||A*x|| ≈ |λ|
        lambda = y_norm;

        // Normalize y for next iteration.
        for (int i = 0; i < n; i++){
            x[i] = y[i] / y_norm;
        }
    }
    
    free(y);
    return lambda;
}



// double powerIteration_old(const csr_matrix *A,
//                       double *x,        // user-allocated vector, size A->n
//                       int max_iter,
//                       double tol)
// {
//     int n = A->n;
//     double *y = (double *)malloc(n * sizeof(double));
//     if(!y) {
//         fprintf(stderr, "Allocation error in powerIteration.\n");
//         exit(1);
//     }

//     // We start with a random or user-provided x. (Here, let's do random.)
//     srand((unsigned)time(NULL));
//     for(int i = 0; i < n; i++) {
//         x[i] = (double)rand() / RAND_MAX;
//     }

//     // We'll track the old and new eigenvalue estimates.
//     double lambda_old = 1.0;
//     double lambda_new = 0.0;

//     // Start iteration
//     for(int iter = 1; iter <= max_iter; iter++) {
//         // (a) y = A * x
//         matVec(A, x, y);

//         // (b) find largest abs value in y -> that's our lambda_new
//         lambda_new = fabs(y[0]);
//         for(int i = 1; i < n; i++) {
//             double val_abs = fabs(y[i]);
//             if(val_abs > lambda_new) {
//                 lambda_new = val_abs;
//             }
//         }

//         // (c) normalize y by dividing each entry by lambda_new
//         //     to get the next x
//         if(lambda_new < 1e-14) {
//             // x is nearly the zero vector => break
//             fprintf(stderr, "Warning: x is near zero. Stopping.\n");
//             break;
//         }
//         for(int i = 0; i < n; i++) {
//             x[i] = y[i] / lambda_new;
//         }

//         // (d) check for convergence
//         if(fabs(lambda_new - lambda_old) < tol) {
//             printf("Converged at iteration %d\n", iter);
//             break;
//         }

//         // (e) update lambda_old
//         lambda_old = lambda_new;
//     }

//     free(y);
//     return lambda_new;
// }


// Power iteration that returns the eigenvector and prints the largest eigenvalue
double* PIV(csr_matrix* A, double *x, int max_iter) {
    int n = A->n;
    double *y = (double*)malloc(n * sizeof(double));
    if (y == NULL) {
        fprintf(stderr, "Memory allocation error.\n");
        exit(EXIT_FAILURE);
    }

    // Randomize x
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        x[i] = (double)rand() / (double)RAND_MAX;
    }

    // Normalize initial x
    {
        double norm_x = vector_norm(n, x);
        if (norm_x > 0) {
            for (int i = 0; i < n; i++) {
                x[i] /= norm_x;
            }
        }
    }

    // Power iteration loop
    for (int iter = 0; iter < max_iter; iter++) {
        // y = A*x
        matVec(A, x, y);

        // Compute norm of y and normalize
        double y_norm = vector_norm(n, y);
        if (y_norm > 1e-6) {
            for (int i = 0; i < n; i++) {
                x[i] = y[i] / y_norm;
            }
        } else {
            // Handle a zero vector or near-zero norm
            fprintf(stderr, "Warning: x is near zero. Stopping.\n");
            break;
        }
    }

    // After final iteration, approximate largest eigenvalue via Rayleigh quotient
    // lambda ~ x^T A x
    matVec(A, x, y); 
    double lambda = dot_product(n, x, y);
    printf("Approx. Largest Eigenvalue: %f\n", lambda);

    free(y);
    return x;

}







int main(int argc,char ** argv){


    //test the readMMF function
    csr_matrix A;
    //readMMF("data/cfd1/cfd1.mtx",&A);
     readMMF("data/bcsstk03/bcsstk03.mtx",&A);
    //readMMF("test.mtx",&A);
    printf("Matrix size: %d x %d\n",A.n,A.m);
    //let test the power iteration method
    double x[A.n];
   // readVec("vector.txt",x,A.n);
    int max_iter = 50;
    //double lambda = PowerIteration(&A,x,max_iter);
    //double lambda = powerIteration_simplified(&A,x,max_iter);
    //double lambda = powerIteration_old(&A,x,max_iter,1e-6);
    
    
  //  printf("Largest eigenvalue: %lf\n",lambda);
    printf("--------------------------------\n");
    double *eigenvector = PIV(&A, x, max_iter);
    printf("Eigen vector found\n"); 
   




   
   



    return 0;
}