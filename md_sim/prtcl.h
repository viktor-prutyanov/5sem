#ifndef PRTCL_H
#define PRTCL_H

#include "vector3.h"

struct prtcl_t
{
    struct vector3_t r;
    struct vector3_t v;
    struct vector3_t r_eff;
};

struct prtcl_t *prtcl_resolve_eff_r(struct prtcl_t *p, double size);
double prtcl_K(struct prtcl_t *p);
struct vector3_t *prtcl_r(struct vector3_t *r, struct prtcl_t *p1, struct prtcl_t *p2);
struct vector3_t *prtcl_r_eff(struct vector3_t *r, struct prtcl_t *p1, struct prtcl_t *p2, double size);
double prtcl_r2(struct prtcl_t *p1, struct prtcl_t *p2);
double prtcl_r2_eff(struct prtcl_t *p1, struct prtcl_t *p2, double size);

#endif /* PRTCL_H */
