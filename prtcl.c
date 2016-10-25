#include "prtcl.h"
#include <math.h>

struct prtcl_t *prtcl_resolve_eff_r(struct prtcl_t *p, double size)
{
#define EFF(i) \
    if (p->r.i > 0)                                             \
        p->r_eff.i = fmod(p->r.i + size / 2, size) - size / 2;  \
    else                                                        \
        p->r_eff.i = fmod(p->r.i - size / 2, size) + size / 2;
    EFF(x);
    EFF(y);
    EFF(z);
#undef EFF

    return p;
}

double prtcl_K(struct prtcl_t *p)
{
    return vector3_square(&p->v) / 2;
}

struct vector3_t *prtcl_r(struct vector3_t *r, struct prtcl_t *p1, struct prtcl_t *p2)
{
#define DELTA(i) r->i = p1->r.i - p2->r.i;
    DELTA(x);
    DELTA(y);
    DELTA(z);
#undef DELTA
    return r;
}

struct vector3_t *prtcl_r_eff(struct vector3_t *r, struct prtcl_t *p1, struct prtcl_t *p2, double size)
{
#define DELTA(i) \
    r->i = p1->r_eff.i - p2->r_eff.i;           \
    if (r->i > size / 2)                        \
        r->i -= size;                           \
    else if (r->i < -size / 2)                  \
        r->i += size;
    DELTA(x);
    DELTA(y);
    DELTA(z);
#undef DELTA
    return r;
}

double prtcl_r2(struct prtcl_t *p1, struct prtcl_t *p2)
{
    struct vector3_t r;
    prtcl_r(&r, p1, p2);
    return vector3_square(&r);
}

double prtcl_r2_eff(struct prtcl_t *p1, struct prtcl_t *p2, double size)
{
    struct vector3_t r;
    prtcl_r_eff(&r, p1, p2, size);
    return vector3_square(&r);
}
