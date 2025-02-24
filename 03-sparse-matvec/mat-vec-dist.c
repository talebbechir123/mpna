#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mmio.h"
#include <math.h>
#include <mpi.h>




int main(int argc, char**argv){
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
   // printf("Hello from %d of %d\n", rank, size);
   printf("Hello from %d of %d\n", rank, size);
    
    //finish MPI
    MPI_Finalize();
    

    return 0;
}

