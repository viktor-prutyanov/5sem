#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "MDSubDomain.h"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int mpiSize, mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    double dt = 0.001;
    double density = 1.0;
    double targetT = 1.0;
    double time = 1;
    unsigned long int N = atol(argv[1]);
    double L = 16;//pow(N / density, 1./3);

    if (mpiRank == 0)
        printf("size = %lf, N = %ld, time = %lg, dt = %lg, target_temp = %lg, density = %lg\n", 
            L, N, time, dt, targetT, density);

    double t = 0.;

    MDSubDomain mdSubDomain(N, L, mpiSize, mpiRank, targetT);
    while (t < time)
    {
        mdSubDomain.Step(dt, t < time / 2);
        t += dt;
        if (mpiRank == 0)
        {
            double K = mdSubDomain.K;
            double U = mdSubDomain.U;
            printf("t = %lg, <K> = %lg, <U> = %lg, <E> = %lg\n", t, K / N, U / N, (K + U) / N);
        }
    }

    MPI_Finalize();
    return 0;
}

