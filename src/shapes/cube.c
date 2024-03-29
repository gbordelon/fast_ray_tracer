#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "shapes.h"
#include "bounding_box.h"

#include "cube.h"

struct two_doubles {
    double a[2];
};

struct two_doubles
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
        tmin = tmin_numerator * INFINITY;
        if (isnan(tmin)) {
            tmin = INFINITY;
            if (tmin_numerator < 0) {
                tmin = -INFINITY;
            }
        }
        tmax = tmax_numerator * INFINITY;
        if (isnan(tmax)) {
            tmax = INFINITY;
            if (tmax_numerator < 0) {
                tmax = -INFINITY;
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
cube_local_intersect(Shape cube, Ray r, bool stop_after_first_hit)
{
    struct two_doubles xtmin_xtmax = check_axis(r->origin[0], r->direction[0]);
    struct two_doubles ytmin_ytmax = check_axis(r->origin[1], r->direction[1]);
    struct two_doubles ztmin_ztmax = check_axis(r->origin[2], r->direction[2]);

    double tmin = fmax(fmax(xtmin_xtmax.a[0], ytmin_ytmax.a[0]), ztmin_ztmax.a[0]);
    double tmax = fmin(fmin(xtmin_xtmax.a[1], ytmin_ytmax.a[1]), ztmin_ztmax.a[1]);

    if (tmin > tmax) {
        return NULL;
    }

    Intersections xs = cube->xs;
    xs->num = 0;
    Intersection x = xs->xs;
    intersection(tmin, cube, x++);
    intersection(tmax, cube, x);
    xs->num = 2;

    return xs;
}

void
cube_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    double abs_x = fabs(local_point[0]);
    double abs_y = fabs(local_point[1]);
    double abs_z = fabs(local_point[2]);
    double maxc = fmax(fmax(abs_x, abs_y), abs_z);

    vector_default(res);

    if (equal(maxc, abs_x)) {
        res[0] = local_point[0];
    } else if (equal(maxc, abs_y)) {
        res[1] = local_point[1];
    } else {
        res[2] = local_point[2];
    }
}

void
cube(Shape s)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_CUBE;
    bounding_box(&(s->bbox));
    s->bbox_valid = false;
    s->xs = intersections_empty(2);

    s->intersect = shape_intersect;
    s->local_intersect = cube_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = cube_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = shape_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}

Shape
cube_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    cube(s);
    return s;
}
