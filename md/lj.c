#include "lj.h"
#include "prtcl.h"

#define RC_2 16
#define U0_3 0.00547944174
#define U0_4 0.00097632

double lj_U(struct prtcl_t *p1, struct prtcl_t *p2, double size)
{
    double r2 = prtcl_r2_eff(p1, p2, size);
    if (r2 > RC_2)
        return 0;
    double r6 = r2 * r2 * r2;
    return 4 * (1 / r6 - 1) / r6 + U0_4;
}

struct vector3_t *lj_F(struct vector3_t *F, struct prtcl_t *p1, struct prtcl_t *p2, double size)
{
    double r2 = prtcl_r2_eff(p1, p2, size);
    if (r2 > RC_2)
        return vector3_zero(F);
    double r6 = r2 * r2 * r2;
    prtcl_r_eff(F, p1, p2, size);
    return vector3_mul(F, 24 * (2 / r6 - 1) / (r6 * r2));
}
