#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "cube.h"
#include "shapes.h"

struct two_doubles {
    double a[2];
};

struct two_doubles // TODO does this work? 
check_axis(double origin, double direction)
{
    struct two_doubles min_max;

    double tmin_numerator = -1 - origin;
    double tmax_numerator = 1 - origin;
    double tmin, tmax;

    if (fabs(direction) >= EPSILON) {
        tmin = tmin_numerator / direction;
        tmax = tmax_numerator / direction;
    } else {
        tmin = tmin_numerator * DBL_MAX;
        if (isnan(tmin)) {
            tmin = DBL_MAX;
            if (tmin_numerator < 0) {
                tmin = -tmin;
            }
        }
        tmax = tmax_numerator * DBL_MAX;
        if (isnan(tmax)) {
            if (tmax_numerator < 0) {
                tmax = -tmax;
            }
        }
    }

    if (tmin > tmax) {
        min_max.a[0] = tmax;
        min_max.a[1] = tmin;
    } else {
        min_max.a[0] = tmin;
        min_max.a[1] = tmax;
    }

    return min_max;
}

Intersections
cube_local_intersect(Shape cube, Ray r)
{
    struct two_doubles xtmin_xtmax = check_axis(r->origin[0], r->direction[0]);
    struct two_doubles ytmin_ytmax = check_axis(r->origin[1], r->direction[1]);
    struct two_doubles ztmin_ztmax = check_axis(r->origin[2], r->direction[2]);

    double tmin = fmax(fmax(xtmin_xtmax.a[0], ytmin_ytmax.a[0]), ztmin_ztmax.a[0]);
    double tmax = fmin(fmin(xtmin_xtmax.a[1], ytmin_ytmax.a[1]), ztmin_ztmax.a[1]);

    if (tmin > tmax) {
        return intersections_empty(0);
    }

    Intersections xs = intersections_empty(2);
    Intersection x = xs->xs;
    intersection(tmin, cube, x++);
    intersection(tmax, cube, x);
    xs->num = 2;

    return xs;
}

Vector
cube_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    double abs_x = fabs(local_point->arr[0]);
    double abs_y = fabs(local_point->arr[1]);
    double abs_z = fabs(local_point->arr[2]);
    double maxc = fmax(fmax(abs_x, abs_y), abs_z);

    if (equal(maxc, abs_x)) {
        return vector(local_point->arr[0], 0, 0);
    } else if (equal(maxc, abs_y)) {
        return vector(0, local_point->arr[1], 0);
    }

    return vector(0,0,local_point->arr[2]);
}

void
cube(Shape s)
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_CUBE;

    s->intersect = shape_intersect;
    s->local_intersect = cube_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = cube_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;
}

Shape
cube_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    cube(s);
    return s;
}
