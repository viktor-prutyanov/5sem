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
    unsigned long int N = atol(argv[1]);
    double L = mpiSize * Rc;
    double density = N / pow(L, 3);
    double targetT = 1.0;
    double time = 300 * dt;
    FILE *dump_file; 
        
    if (mpiRank == 0)
        dump_file = fopen("dump.csv", "w");

    if (mpiRank == 0)
        printf("L = %lg, N = %ld, time = %lg, dt = %lg, targetT = %lg, density = %lg\n", 
            L, N, time, dt, targetT, density);

    double t = 0.;

    MDSubDomain sd(N, L, mpiSize, mpiRank, targetT);
    //sd.MicroDump(stdout);

    while (t < time)
    {
        sd.Step(dt, t < time / 2);
        t += dt;
        if (mpiRank == 0)
            sd.MacroDump(dump_file, t);
    }

    if (mpiRank == 0)
    {
        printf("t = %lg, <K> = %lg, <U> = %lg, <E> = %lg\n", t, sd.K / N, sd.U / N, (sd.K + sd.U) / N); 
        fclose(dump_file);
    }

    MPI_Finalize();
    return 0;
}

