#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank );
    double data = rank * rank;
    printf("Process %d / %d data = %lg (initial)\n", rank, size, data);
    
    MPI_Status status;
    if (rank % 2 == 0)
    {
        MPI_Send(&data, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&data, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
    }
    else
    {
        MPI_Send(&data, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&data, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
    }

    printf("Process %d / %d data = %lg (final)\n", rank, size, data);
    
    MPI_Finalize();
    return 0;
}
