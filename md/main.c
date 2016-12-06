#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "md_sim.h"
#include "get_type_opt.h"

#define CRYSTAL_RND 0.05

int main(int argc, char *argv[])
{
    double dt = 0.;
    double density = 0.;
    double target_T = 0.;
    double size = 0.;
    double time = 0.;
    double dv = 0.0001;
    long int N = 0;
    FILE *dump_file = NULL;

    int opt = 0;
    
    const char usage_str[] = "usage: %s [-s size] [-N particles] [-r density]\
[-t time] [-d delta_t] [-T target_temp] [-o dump_file] [-v dv in distribution]\n";
    
    while ((opt = getopt(argc, argv, "s:N:t:d:o:T:r:v:")) != -1)
    {
        switch(opt)
        {
        case 's':
            size = get_pos_double_opt(optarg, "size of md_sim");
            break;
        case 'N':
            N = get_pos_long_opt(optarg, "particles number");
            break;
        case 't':
            time = get_pos_double_opt(optarg, "time to simulate");
            break;
        case 'd':
            dt = get_pos_double_opt(optarg, "time delta");
            break;
        case 'T':
            target_T = get_pos_double_opt(optarg, "target temperature");
            break;
        case 'v':
            dv = get_pos_double_opt(optarg, "dv in distribution");
            break;
        case 'r':
            density = get_pos_double_opt(optarg, "density");
            break;
        case 'o':
            dump_file = fopen(optarg, "w");
            if (!dump_file)
            {
                fprintf(stderr, "Args error: dump file opening failed.\n");
                return -1;
            }
            break;
        default:
            fprintf(stderr, usage_str, argv[0]);
            return -1;
            break;
        }
    }
    
    if (optind != argc)
    {
        fprintf(stderr, usage_str, argv[0]);
        return -1;
    }


    if ((density == 0.) && (N != 0) && (size != 0))
        density = N / (size*size*size);
    else if ((size == 0.) && (density != 0) && (N != 0))
    {
        size = pow(N / density, 1./3);
    }
    else if ((N == 0) && (size != 0) && (density != 0.))
        N = (long int)round(density * size*size*size);
    else
    {
        fprintf(stderr, "Args error: enter as args strictly 2 values from list: density, N, size.\n");
        return -1;
    }

    printf("size = %lf, N = %ld, time = %lg, dt = %lg, target_temp = %lg, density = %lg, dv = %lg\n", 
        size, N, time, dt, target_T, density, dv);

    struct md_sim_t md_sim;
    md_sim_ctor(&md_sim, N, size, dt, target_T, dump_file);
    md_sim_init_crystal(&md_sim, CRYSTAL_RND);
    md_sim_run(&md_sim, time);
    md_sim_dump_v_distribution(&md_sim, dv, "distr.csv");
    md_sim_dtor(&md_sim);

    fclose(dump_file);

    return 0;
}

