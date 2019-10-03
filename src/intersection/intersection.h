#ifndef INTERSECTION
#define INTERSECTION

#include <stdbool.h>

struct shape;
typedef struct shape *Shape;

typedef struct intersection {
    double t;
    double u;
    double v;
    Shape object;
} *Intersection;

typedef struct intersections {
    Intersection xs;
    size_t array_len;
    size_t num;
} *Intersections;

void intersection(const double t, const Shape sh, Intersection x);
Intersection intersection_alloc(const double t, const Shape sh);

void intersection_with_uv(const double t, const double u, const double v, const Shape sh, Intersection x);
Intersection intersection_with_uv_alloc(const double t, const double u, const double v, const Shape sh);

Intersection hit(const Intersections xs, const bool filter_shadow_casters);
Intersections intersections_empty(const size_t num);

void intersections_free(Intersections xs);

void intersections_sort(Intersections xs);
void intersections_reverse(Intersections xs);

int intersection_to_string(char *buf, size_t n, Intersection x);


#endif
