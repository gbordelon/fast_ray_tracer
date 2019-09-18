#include <math.h>
#include <stdlib.h>

#include "sphere.h"
#include "shapes.h"

Intersections
sphere_local_intersect(Shape sphere, Ray r)
{
    double sphere_origin[4] = { 0.0, 0.0, 0.0, 0.0 };;
    Vector sphere_to_ray = vector_from_arrays_alloc(r->origin, sphere_origin);
    double a = array_dot(r->direction, r->direction);
    double b = 2 * array_dot(r->direction, sphere_to_ray->arr);
    double c = vector_dot(sphere_to_ray, sphere_to_ray) - 1.0;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return intersections_empty(0);
    }
    discriminant = sqrt(discriminant);
    a = 1.0 / (2 * a);

    Intersections xs = intersections_empty(2);
    Intersection x = xs->xs;
    intersection((-b - discriminant) * a, sphere, x++);
    intersection((-b + discriminant) * a, sphere, x);
    xs->num = 2;

    vector_free(sphere_to_ray);

    return xs;
}

Vector
sphere_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    Point origin = point(0,0,0);
    Vector v = vector_from_points_alloc(local_point, origin);

    point_free(origin);

    return v;
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

    s->intersect = shape_intersect;
    s->local_intersect = sphere_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = sphere_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;
}

Shape
sphere_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    sphere(s);
    return s;
}
