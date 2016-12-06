#include <string.h>
#include <strings.h>
#include <stdio.h>

#include "vector3.h"

struct vector3_t *vector3_print(const char *prefix, struct vector3_t *v, const char *suffix)
{
    printf(prefix);
    printf("{%lg; %lg; %lg}", v->x, v->y, v->z);
    printf(suffix);
    return v;
}

struct vector3_t *vector3_zero(struct vector3_t *v)
{
    bzero(v, sizeof(struct vector3_t));
    return v;
}

double vector3_square(struct vector3_t *v)
{
    return v->x * v->x + v->y * v->y + v->z * v->z;
}

struct vector3_t *vector3_mul(struct vector3_t *v, double k)
{
    v->x = v->x * k;
    v->y = v->y * k;
    v->z = v->z * k;
    return v;
}

struct vector3_t *vector3_add(struct vector3_t *v, struct vector3_t *v1, struct vector3_t *v2)
{
    v->x = v1->x + v2->x;
    v->y = v1->y + v2->y;
    v->z = v1->z + v2->z;
    return v;
}

struct vector3_t *vector3_copy(struct vector3_t *v_dst, struct vector3_t *v_src)
{
    return memcpy(v_dst, v_src, sizeof(struct vector3_t));
}
