#include <math.h>
#include <stdlib.h>

#include "toroid.h"
#include "shapes.h"

/*
Intersections
toroid_local_intersect(Shape toroid, Ray r)
{
    double toroid_origin[4] = { 0.0, 0.0, 0.0, 0.0 };
    // does this makes sense for a toroid?
    Vector toroid_to_ray = vector_from_arrays_alloc(r->origin, toroid_origin);

    double a = array_dot(r->direction, r->direction);
    double b = 2 * array_dot(r->direction, toroid_to_ray->arr);
    double c = vector_dot(toroid_to_ray, toroid_to_ray) - 1.0;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return intersections_empty(0);
    }
    discriminant = sqrt(discriminant);
    a = 1.0 / (2 * a);

    Intersections xs = intersections_empty(2);
    Intersection x = xs->xs;
    intersection((-b - discriminant) * a, toroid, x++);
    intersection((-b + discriminant) * a, toroid, x);
    xs->num = 1;

    vector_free(toroid_to_ray);

    return xs;
}

Vector
toroid_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    double toroid_origin[4] = {0.0, 0.0, 0.0, 0.0};

    Vector v1 = vector_from_arrays_alloc(local_point->arr, toroid_origin);
    Vector v2 = vector_normalize_alloc(v1);
    vector_scale(v2, sh->fields.toroid.radius1);


    vector_from_arrays(local_point->arr, v2->arr, v1); // hack
    vector_normalize(v1, v2);

    vector_free(v1);

    return v2;
}

void
toroid(Shape s)
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_TOROID;
    s->fields.toroid.radius1 = 0.5;
    s->fields.toroid.radius2 = 0.5;

    s->intersect = shape_intersect;
    s->local_intersect = toroid_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = toroid_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;
}

Shape
toroid_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    toroid(s);
    return s;
}

*/
