#ifndef LJ_H
#define LJ_H

#include "prtcl.h"

double lj_U(struct prtcl_t *p1, struct prtcl_t *p2, double size);
struct vector3_t *lj_F(struct vector3_t *F, struct prtcl_t *p1, struct prtcl_t *p2, double size);

#endif /* LJ_H */
