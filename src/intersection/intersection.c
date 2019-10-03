#include <stdbool.h>

#include "../shapes/shapes.h"

#include "intersection.h"

void
intersection(const double t, const Shape sh, Intersection x)
{
    x->t = t;
    x->object = sh;
    x->u = -1;
    x->v = -1;
}

void
intersection_with_uv(const double t, const double u, const double v, const Shape sh, Intersection x)
{
    x->t = t;
    x->object = sh;
    x->u = u;
    x->v = v;
}

Intersection
intersection_alloc(const double t, const Shape sh)
{
    Intersection x = (Intersection) malloc(sizeof(struct intersection));
    intersection(t, sh, x);
    return x;
}

Intersection
intersection_uv_alloc(const double t, const double u, const double v, const Shape sh)
{
    Intersection x = (Intersection) malloc(sizeof(struct intersection));
    intersection_with_uv(t, u, v, sh, x);
    return x;
}

Intersection
hit(const Intersections xs, const bool filter_shadow_casters)
{
    int i;
    Intersection x;

    for (i = 0, x = xs->xs; i < xs->num; i++, x++) {
        if (x->t > 0 && (!filter_shadow_casters || (filter_shadow_casters && x->object->material->casts_shadow))) {
            return x;
        }
    }

    return NULL;
}

Intersections
intersections_empty(const size_t num)
{
    Intersections xs = (Intersections) malloc(sizeof(struct intersections));
    if (num > 0) {
        xs->xs = (Intersection) malloc(num * sizeof(struct intersection));
    } else {
        xs->xs = NULL;
    }

    xs->array_len = num;
    xs->num = 0;

    return xs;
}

void
intersection_free(Intersection x)
{
    if (x != NULL) {
        free(x);
    }
}

void
intersections_free(Intersections xs)
{
    if (xs != NULL) {
        if (xs->xs != NULL) {
            free(xs->xs);
            xs->xs = NULL;
        }
        free(xs);
    }
}

int
intersection_to_string(char *buf, size_t n, Intersection x)
{
    int used = snprintf(buf, n, "Intersection:\n\tt: %f\n\tu: %f\n\tv: %f\n",
                        x->t,
                        x->u,
                        x->v);
    return used + shape_to_string(buf + used, n - used, x->object);
}

int
sort_intersections_asc(const void *p, const void *q)
{
    double l = ((Intersection)p)->t;
    double r = ((Intersection)q)->t;
    if (l - r < 0) {
        return -1;
    } else if (l - r > 0) {
        return 1;
    }
    return 0;
}

int
sort_intersections_desc(const void *p, const void *q)
{
    double l = ((Intersection)p)->t;
    double r = ((Intersection)q)->t;
    if (l - r < 0) {
        return 1;
    } else if (l - r > 0) {
        return -1;
    }
    return 0;
}

void
intersections_reverse(Intersections xs)
{
    // sort xs by xs->xs->t descending
    qsort((void*)xs->xs, xs->num, sizeof(struct intersection), sort_intersections_desc);
}

void
intersections_sort(Intersections xs)
{
    // sort xs by xs->xs->t ascending
    qsort((void*)xs->xs, xs->num, sizeof(struct intersection), sort_intersections_asc);
}
