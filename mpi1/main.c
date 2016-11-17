#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank );
    double data = (rank == 0) ? 1.234 : 0;
    //printf("Process %d / %d data = %lg (initial)\n", rank, size, data);
    
    if (rank == 0)
    {
        MPI_Send(&data, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    }
    else if (rank == size - 1)
    {
        MPI_Status status;
        MPI_Recv(&data, 1, MPI_DOUBLE, size - 2, 0, MPI_COMM_WORLD, &status);
    }
    else
    {
        MPI_Status status;
        MPI_Recv(&data, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&data, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    }
        
    printf("Process %d / %d data = %lg (final)\n", rank, size, data);
    
    MPI_Finalize();
    return 0;
}
