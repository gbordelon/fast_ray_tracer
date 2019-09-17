#include <stdlib.h>
#include <string.h>

#include "shapes.h"

Intersection
hit(Intersections xs)
{
    int i;
    Intersection x;

    for (i = 0, x = xs->xs; i < xs->len; i++, x++) {
        if (x->t > 0) {
            return x;
        }
    }

    return NULL;
}

Intersections
intersections_empty(size_t num)
{
    Intersections xs = (Intersections) malloc(sizeof(struct intersections));
    xs->xs = (Intersection) malloc(num * sizeof(struct intersection));
    xs->len = num;

    return xs;
}

void
intersections_free(Intersections xs)
{
    if (xs != NULL) {
        if (xs->xs != NULL) {
            free(xs);
        }
        free(xs);
    }
}


Ray
ray_alloc(Point origin, Vector direction)
{
    Ray ray = (Ray) malloc(sizeof(struct ray));
    memcpy(ray->origin, origin->arr, sizeof(ray->origin));
    memcpy(ray->direction, direction->arr, sizeof(ray->direction));

    return ray;
}

Material
material()
{
    Material m = (Material) malloc(sizeof(struct material));
    m->color[0] = 1.0;
    m->color[1] = 1.0;
    m->color[2] = 1.0;
    m->ambient = 0.1;
    m->diffuse = 0.9;
    m->specular = 0.9;
    m->shininess = 200.0;
    m->reflective = 0.0;
    m->transparency = 0.0;
    m->refractive_index = 1.0;
    m->casts_shadow = true;
    m->pattern = NULL;

    return m;
}

Intersections
shape_intersect(Shape sh, Ray r);

Vector
shape_normal_at(Shape sh, Point world_point, Intersection hit);

Intersections
sphere_local_intersect(Shape sphere, Ray r);

Vector
sphere_local_normal_at(Shape sh, Point local_point, Intersection hit);

Vector
shape_normal_to_world(Shape sh, Vector local_normal);

Point
shape_world_to_object(Point pt);

void
shape_divide(Shape sh, size_t threshold);

bool
shape_includes(Shape a, Shape b);

Shape
sphere()
{
    Shape s = (Shape) malloc(sizeof(struct shape));

    s->transform = matrix_identity_alloc();
    s->material = material();
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

    return s;
}
