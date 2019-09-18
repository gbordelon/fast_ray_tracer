#include <math.h>
#include <stdlib.h>

#include "linalg.h"
#include "plane.h"
#include "shapes.h"

Intersections
plane_local_intersect(Shape plane, Ray r)
{
    if (fabs(r->direction[1]) < EPSILON) {
        return intersections_empty(0);
    }

    Intersections xs = intersections_empty(1);
    Intersection x = xs->xs;
    intersection(-r->origin[1] / r->direction[1], plane, x);
    xs->num = 1;

    return xs;
}

Vector
plane_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    return vector(0,1,0);
}

void
plane(Shape s)
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_PLANE;

    s->intersect = shape_intersect;
    s->local_intersect = plane_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = plane_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;
}

Shape
plane_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    plane(s);
    return s;
}
