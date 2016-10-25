#ifndef VECTOR3_H
#define VECTOR3_H

struct vector3_t
{
    double x;
    double y;
    double z;
};

double vector3_square(struct vector3_t *v);
struct vector3_t *vector3_print(const char *prefix, struct vector3_t *v, const char *suffix);
struct vector3_t *vector3_zero(struct vector3_t *v);
struct vector3_t *vector3_mul(struct vector3_t *v, double k);
struct vector3_t *vector3_add(struct vector3_t *v, struct vector3_t *v1, struct vector3_t *v2);
struct vector3_t *vector3_copy(struct vector3_t *v_dst, struct vector3_t *v_src);

#endif /* VECTOR3_H */
