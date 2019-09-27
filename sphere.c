#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "sphere.h"
#include "shapes.h"

Intersections
sphere_local_intersect(Shape sphere, Ray r)
{
    double sphere_origin[4] = { 0.0, 0.0, 0.0, 0.0 };
    struct v sphere_to_ray;
    vector_from_arrays(r->origin, sphere_origin, &sphere_to_ray);
    double a = array_dot(r->direction, r->direction);
    double b = 2 * array_dot(r->direction, sphere_to_ray.arr);
    double c = vector_dot(&sphere_to_ray, &sphere_to_ray) - 1.0;
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
    memcpy(res->arr, local_point->arr, 3 * sizeof(double));
    res->arr[3] = 0.0;
}

void
sphere(Shape s)
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_SPHERE;
    s->bbox = NULL;
    s->bbox_inverse = NULL;
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
