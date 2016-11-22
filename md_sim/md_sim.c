#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>

#include "vector3.h"
#include "md_sim.h"
#include "prtcl.h"
#include "lj.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define K2T(K,N) (K * 2./ (N * 3.))
#define TAU 0.1

int md_sim_ctor(struct md_sim_t *md_sim, long int n_prtcls, double size, double dt, double target_T, FILE *dump_file)
{
    md_sim->prtcls = (struct prtcl_t *)calloc(n_prtcls, sizeof(struct prtcl_t));
    if (!md_sim->prtcls)
        return -1;
    
    md_sim->size = size;
    md_sim->n = n_prtcls;
    md_sim->dt = dt;
    md_sim->dump_file = dump_file;

    md_sim->targets.K = n_prtcls * target_T * 3./2;
    
    md_sim->tweaks.tau = TAU;
    md_sim->tweaks.thermostat = 1;

    return 0;
}

int md_sim_dtor(struct md_sim_t *md_sim)
{
    md_sim->n = 0;
    free(md_sim->prtcls);
    return 0;
}

int md_sim_init_crystal(struct md_sim_t *md_sim, double rnd_param)
{
    int q = (int)round(pow(md_sim->n, 1./3));
    double step = md_sim->size / q;
    
    srand(time(NULL));
    
    long int count = 0;
    for (int i = 0; i < q; ++i)
        for (int j = 0; j < q; ++j)
            for (int k = 0; k < q; ++k)
            {
                md_sim->prtcls[count].r.x = step / 2 + i * step + step * rnd_param * (rand() * 1. - RAND_MAX/2) / RAND_MAX;
                md_sim->prtcls[count].r.y = step / 2 + j * step + step * rnd_param * (rand() * 1. - RAND_MAX/2) / RAND_MAX;
                md_sim->prtcls[count].r.z = step / 2 + k * step + step * rnd_param * (rand() * 1. - RAND_MAX/2) / RAND_MAX;
                ++count;
            }
    
    for (long int i = count; i < md_sim->n; ++i)
    {
        md_sim->prtcls[i].r.x = rand() * md_sim->size * 1. / RAND_MAX;
        md_sim->prtcls[i].r.y = rand() * md_sim->size * 1. / RAND_MAX;
        md_sim->prtcls[i].r.z = rand() * md_sim->size * 1. / RAND_MAX;
    }
    
    for (long int i = 0; i < md_sim->n; ++i)
        memset(&md_sim->prtcls[i].v, 0, 3 * sizeof(double));
    
    md_sim_print_prtcls(md_sim);
    fprintf(md_sim->dump_file, "t,U,K,E,<U>,<K>,<E>,T,lambda,tau,t_on\n");
    
    md_sim_resolve_eff_r(md_sim);
    md_sim_update_stats(md_sim);
    md_sim_dump_stats(md_sim, -md_sim->dt);       

    return 0;
}

int md_sim_dump_v_distribution(struct md_sim_t *md_sim, double dv, const char *filename)
{
    FILE *file = fopen(filename, "w");

    long int *v_abs_norm = (long int *)calloc(md_sim->n, sizeof(long int));

    long int v_max = 0;
    for (long int i = 0; i < md_sim->n; ++i)
    {
        v_abs_norm[i] = lrint(sqrt(vector3_square(&md_sim->prtcls[i].v)) / dv);
        if (v_abs_norm[i] > v_max) 
            v_max = v_abs_norm[i];
    }

    fprintf(file, "V_abs,dN,V_maxwell\n");
    
    for (long int v = 0; v < v_max; ++v)
    {
        long int count = 0;
        for (long int i = 0; i < md_sim->n; ++i)
            if (v_abs_norm[i] == v)
                ++count;
        fprintf(file, "%ld,%lg,%lg\n", v, count * 1./ md_sim->n, 
                maxwell_by_abs(1, 1, v * dv, K2T(md_sim->targets.K, md_sim->n)) * dv);
    }

    free(v_abs_norm);
    return 0;
}

double maxwell_by_abs(double k, double m, double v, double T)
{
   return 4 * M_PI * v*v * pow(m / (2 * M_PI * k * T), 1.5) * exp(- m * v*v / (2 * k * T));
}

int md_sim_dump_stats(struct md_sim_t *md_sim, double t)
{
    return fprintf(md_sim->dump_file, 
        "%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%d\n", t, 
        md_sim->stats.U, 
        md_sim->stats.K, 
        md_sim->stats.U + md_sim->stats.K, 
        md_sim->stats.U / md_sim->n, 
        md_sim->stats.K / md_sim->n, 
        (md_sim->stats.U + md_sim->stats.K) / md_sim->n,
        K2T(md_sim->stats.K, md_sim->n),
        md_sim->stats.lambda,
        md_sim->tweaks.tau,
        md_sim->tweaks.thermostat);
}

void md_sim_print_prtcls(struct md_sim_t *md_sim)
{
    printf("{\n");
    for (long int i = 0; i < md_sim->n; ++i)
    {
        printf("\t[%ld]: ", i);
        vector3_print("r = ", &md_sim->prtcls[i].r, " ");
        vector3_print("v = ", &md_sim->prtcls[i].v, "\n");
    }
    printf("}\n");
}

struct md_sim_stats_t *md_sim_update_stats(struct md_sim_t *md_sim)
{
    md_sim_update_U(md_sim);
    md_sim_update_K(md_sim);
    md_sim_update_lambda(md_sim);
    return &md_sim->stats;
}

void md_sim_lambda_scale(struct md_sim_t *md_sim)
{
    for (long int i = 0; i < md_sim->n; ++i)
        vector3_mul(&md_sim->prtcls[i].v, md_sim->stats.lambda);
}

double md_sim_update_lambda(struct md_sim_t *md_sim)
{
    return (md_sim->stats.lambda = 
        sqrt(1 + md_sim->dt / md_sim->tweaks.tau * (md_sim->targets.K / md_sim->stats.K - 1)));
}

double md_sim_update_U(struct md_sim_t *md_sim)
{
    md_sim->stats.U = 0.;
    for (long int i = 0; i < md_sim->n; ++i)
        for (long int j = i + 1; j < md_sim->n; ++j)
            md_sim->stats.U += lj_U(&md_sim->prtcls[i], &md_sim->prtcls[j], md_sim->size);
    return md_sim->stats.U;
}

double md_sim_update_K(struct md_sim_t *md_sim)
{
    md_sim->stats.K = 0.;
    for (long int i = 0; i < md_sim->n; ++i)
        md_sim->stats.K += prtcl_K(&md_sim->prtcls[i]);
    return md_sim->stats.K;
}

void md_sim_resolve_eff_r(struct md_sim_t *md_sim)
{
    for (long int i = 0; i < md_sim->n ; ++i)
        prtcl_resolve_eff_r(&md_sim->prtcls[i], md_sim->size);
}

int md_sim_run(struct md_sim_t *md_sim, double time_end)
{
    double t = 0.;
    while (t < time_end)
    {
        printf("t = %lg\n", t);
        for (long int i = 0; i < md_sim->n; ++i)
        {
            struct vector3_t dx;
            vector3_copy(&dx, &md_sim->prtcls[i].v);
            vector3_mul(&dx, md_sim->dt);
            vector3_add(&md_sim->prtcls[i].r, &md_sim->prtcls[i].r, &dx);
        }
        
        md_sim_resolve_eff_r(md_sim);

        for (long int i = 0; i < md_sim->n; ++i)
        {
            struct vector3_t F;
            vector3_zero(&F);
                
            for (long int j = 0; j < md_sim->n; ++j)
            {
                if (i == j)
                    continue;
                
                struct vector3_t F_j;
                lj_F(&F_j, &md_sim->prtcls[i], &md_sim->prtcls[j], md_sim->size);

                vector3_add(&F, &F, &F_j);
            }
           
            vector3_mul(&F, md_sim->dt);
            vector3_add(&md_sim->prtcls[i].v, &md_sim->prtcls[i].v, &F);
        }
       
        if (t > time_end/2)
            md_sim->tweaks.thermostat = 0;

        md_sim_update_stats(md_sim);
        md_sim_dump_stats(md_sim, t);       
       
        if (md_sim->tweaks.thermostat)
            md_sim_lambda_scale(md_sim);

        t += md_sim->dt;
    }

    return 0;
}
