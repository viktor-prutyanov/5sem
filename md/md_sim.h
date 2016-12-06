#ifndef MD_SIM_H
#define MD_DIM_H

struct md_sim_targets_t
{
    double K;
};

struct md_sim_stats_t
{
    double K;
    double U;
    double lambda;
};

struct md_sim_tweaks_t
{
    double tau;
    int thermostat;
};

struct md_sim_t
{
    long int n;
    double size;
    struct prtcl_t *prtcls;
    double dt;
    struct md_sim_targets_t targets;
    struct md_sim_stats_t stats;
    struct md_sim_tweaks_t tweaks;
    FILE *dump_file;
};

double maxwell_by_abs(double k, double m, double v, double T);
int md_sim_dump_v_distribution(struct md_sim_t *md_sim, double dv, const char *filename);
void md_sim_lambda_scale(struct md_sim_t *md_sim);
int md_sim_dump_stats(struct md_sim_t *md_sim, double t);
void md_sim_resolve_eff_r(struct md_sim_t *md_sim);
double md_sim_update_lambda(struct md_sim_t *md_sim);
double md_sim_update_K(struct md_sim_t *md_sim);
double md_sim_update_U(struct md_sim_t *md_sim);
struct md_sim_stats_t *md_sim_update_stats(struct md_sim_t *md_sim);
void md_sim_print_prtcls(struct md_sim_t *md_sim);
int md_sim_ctor(struct md_sim_t *md_sim, long int n_prtcls, double size, double dt, double target_T, FILE *dump_file);
int md_sim_dtor(struct md_sim_t *md_sim);
int md_sim_init_crystal(struct md_sim_t *md_sim, double rnd_param);
int md_sim_run(struct md_sim_t *md_sim, double time);

#endif /* MD_SIM_H */
