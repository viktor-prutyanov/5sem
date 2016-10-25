#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "vector3.h"
#include "md_sim.h"
#include "prtcl.h"
#include "lj.h"

#define RND_PARAM 0.05

int md_sim_ctor(struct md_sim_t *md_sim, long int n_prtcls, double size, double dt, double target_T, FILE *dump_file)
{
    md_sim->prtcls = (struct prtcl_t *)calloc(n_prtcls, sizeof(struct prtcl_t));
    if (!md_sim->prtcls)
        return -1;
    
    md_sim->size = size;
    md_sim->n = n_prtcls;
    md_sim->dt = dt;
    md_sim->targets.K = n_prtcls * target_T * 3/2;
    md_sim->dump_file = dump_file;

    return 0;
}

int md_sim_dtor(struct md_sim_t *md_sim)
{
    md_sim->n = 0;
    free(md_sim->prtcls);
    return 0;
}

int md_sim_init(struct md_sim_t *md_sim)
{
    int q = (int)round(pow(md_sim->n, 1./3));
    double step = md_sim->size / q;
    
    srand(time(NULL));
    
    long int count = 0;
    for (int i = 0; i < q; ++i)
        for (int j = 0; j < q; ++j)
            for (int k = 0; k < q; ++k)
            {
                md_sim->prtcls[count].r.x = step / 2 + i * step + step * RND_PARAM * (rand() * 1. - RAND_MAX/2) / RAND_MAX;
                md_sim->prtcls[count].r.y = step / 2 + j * step + step * RND_PARAM * (rand() * 1. - RAND_MAX/2) / RAND_MAX;
                md_sim->prtcls[count].r.z = step / 2 + k * step + step * RND_PARAM * (rand() * 1. - RAND_MAX/2) / RAND_MAX;
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
    fprintf(md_sim->dump_file, "t,U,K,E,<U>,<K>,<E>,lambda\n");
    
    md_sim_resolve_eff_r(md_sim);
    md_sim_update_stats(md_sim);
    md_sim_dump_stats(md_sim, -md_sim->dt);       

    return 0;
}

int md_sim_dump_stats(struct md_sim_t *md_sim, double t)
{
    return fprintf(md_sim->dump_file, "%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e,%015.07e\n", t, 
        md_sim->stats.U, 
        md_sim->stats.K, 
        md_sim->stats.U + md_sim->stats.K, 
        md_sim->stats.U / md_sim->n, 
        md_sim->stats.K / md_sim->n, 
        (md_sim->stats.U + md_sim->stats.K) / md_sim->n, 
        md_sim->stats.lambda);
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

double md_sim_update_lambda(struct md_sim_t *md_sim)
{
    return (md_sim->stats.lambda = 0.);
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

int md_sim_run(struct md_sim_t *md_sim, double time)
{
    double t = 0.;
    while (t < time)
    {
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
        
        md_sim_update_stats(md_sim);
        md_sim_dump_stats(md_sim, t);       

        t += md_sim->dt;
    }

    return 0;
}
