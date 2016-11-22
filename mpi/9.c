#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
    int size, rank;
    if (MPI_Init(&argc, &argv))
        return -1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double data = rank * rank;
    double *recvbuf;
    
    if (rank == 0)
        recvbuf = (double *)malloc(size * sizeof(double));
    
    MPI_Gather(&data, 1, MPI_DOUBLE, recvbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    
    if (rank == 0)
        for(int i = 1; i < size; ++i) 
            printf("Received '%lg' from process %d\n", recvbuf[i], i);                

  MPI_Finalize();
  return 0;
}
