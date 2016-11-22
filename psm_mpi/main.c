/**
 *  Parallel Simpson method
 *
 *  @file main.c
 *  
 *  @date 11.2016
 * 
 *  @copyright GNU GPL v2.0
 * 
 *  @author Viktor Prutyanov mailto:viktor.prutyanov@phystech.edu 
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define FUNC(x) ((x)*exp(x))

#define CORRECT_OR_FINALIZE(st, msg)\
    if (!(st))                      \
    {                               \
        if (rank == 0)              \
            fprintf(stderr, msg);   \
        MPI_Finalize();             \
        return 1;                   \
    }

double calc_segm(double begin, double end, double d, int size, int rank);

int main(int argc, char *argv[])
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    CORRECT_OR_FINALIZE(argc == 4, "Usage: psm [begin] [end] [delta]\n"); 

    double result;
    double subresult = calc_segm(atof(argv[1]), atof(argv[2]), atof(argv[3]), size, rank);
    MPI_Reduce(&subresult, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) 
        printf("Result is %lg\n", result); 

    MPI_Finalize();

    return 0;
}

double calc_segm(double begin, double end, double d, int size, int rank)
{      
    double subsegm_len = (end - begin) / size;
    unsigned long int subsubsegm_num = (unsigned long int)floor(subsegm_len / d);
    double subsegm_begin = begin + rank * subsegm_len;
    double calc_begin = subsegm_begin;
    double result = 0.;

    for (unsigned long int i = 0; i < subsubsegm_num; ++i)
    {
        result += d * (FUNC(calc_begin) + 4 * FUNC(calc_begin + d / 2) + FUNC(calc_begin + d)) / 6;
        calc_begin += d;
    }

    double last_d = subsegm_begin + subsegm_len - calc_begin;
    return result + last_d * (FUNC(calc_begin) + 4 * FUNC(calc_begin + last_d / 2) + FUNC(calc_begin + last_d)) / 6;
}
