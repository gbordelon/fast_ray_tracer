#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../libs/linalg/linalg.h"
#include "../intersection/intersection.h"

#include "shapes.h"
#include "bounding_box.h"

#include "sphere.h"

Intersections
sphere_local_intersect(Shape sphere, Ray r, bool stop_after_first_hit)
{
    Point sphere_origin = { 0.0, 0.0, 0.0, 0.0 };
    Vector sphere_to_ray;
    vector_from_points(r->origin, sphere_origin, sphere_to_ray);
    double a = vector_dot(r->direction, r->direction);
    double b = 2 * vector_dot(r->direction, sphere_to_ray);
    double c = vector_dot(sphere_to_ray, sphere_to_ray) - 1.0;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return NULL;
    }
    discriminant = sqrt(discriminant);
    a = 1.0 / (2 * a);

    Intersections xs = sphere->xs;
    xs->num = 0;

    Intersection x = xs->xs;
    intersection((-b - discriminant) * a, sphere, x++);
    intersection((-b + discriminant) * a, sphere, x);
    xs->num = 2;

    return xs;
}

void
sphere_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    memcpy(res, local_point, sizeof(Vector));
    res[3] = 0.0;
}

void
sphere(Shape s)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_SPHERE;
    bounding_box(&(s->bbox));
    s->bbox_valid = false;
    s->xs = intersections_empty(2);

    s->intersect = shape_intersect;
    s->local_intersect = sphere_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = sphere_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = shape_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}

Shape
sphere_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    sphere(s);
    return s;
}
