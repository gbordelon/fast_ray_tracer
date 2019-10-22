#include <math.h>
#include <stdlib.h>

#include "../libs/linalg/linalg.h"

#include "plane.h"
#include "shapes.h"
#include "bounding_box.h"

Intersections
plane_local_intersect(Shape plane, Ray r, bool stop_after_first_hit)
{
    if (fabs(r->direction[1]) < EPSILON) {
        return NULL;
    }

    Intersections xs = plane->xs;
    xs->num = 0;
    Intersection x = xs->xs;
    intersection(-r->origin[1] / r->direction[1], plane, x);
    xs->num = 1;

    return xs;
}

void
plane_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    res[0] = 0;
    res[1] = 1;
    res[2] = 0;
    res[3] = 0.0;
}

void
plane_bounds(Shape plane, Bounding_box *res)
{
    if (!plane->bbox_valid) {
        plane->bbox_valid = true;
        Point arr = {-INFINITY, 0.0, -INFINITY, 1.0};

        bounding_box_add_array(&(plane->bbox), arr);
        arr[0] = INFINITY;
        arr[1] = 0.0;
        arr[2] = INFINITY;
        bounding_box_add_array(&(plane->bbox), arr);
        bounding_box_transform(&(plane->bbox), plane->transform, &(plane->bbox_inverse));
    }

    *res = plane->bbox;
}

void
plane(Shape s)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_PLANE;
    bounding_box(&(s->bbox));
    s->bbox_valid = false;
    s->xs = intersections_empty(1);

    s->intersect = shape_intersect;
    s->local_intersect = plane_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = plane_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = plane_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}
