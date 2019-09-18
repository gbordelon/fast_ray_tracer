#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shapes.h"

void
intersection(double t, Shape sh, Intersection x)
{
    x->t = t;
    x->object = sh;
    x->u = -1;
    x->v = -1;
}

void
intersection_uv(double t, double u, double v, Shape sh, Intersection x)
{
    x->t = t;
    x->object = sh;
    x->u = u;
    x->v = v;
}

Intersection
intersection_alloc(double t, Shape sh)
{
    Intersection x = (Intersection) malloc(sizeof(struct intersection));
    intersection(t, sh, x);
    return x;
}

Intersection
intersection_uv_alloc(double t, double u, double v, Shape sh)
{
    Intersection x = (Intersection) malloc(sizeof(struct intersection));
    intersection_uv(t, u, v, sh, x);
    return x;
}

Intersection
hit(Intersections xs, bool filter_shadow_casters)
{
    int i;
    Intersection x;

    for (i = 0, x = xs->xs; i < xs->num; i++, x++) {
        if (x->t > 0 && (!filter_shadow_casters || (filter_shadow_casters && x->object->material->casts_shadow))) {
            return x;
        }
    }

    return NULL;
}

Intersections
intersections_empty(size_t num)
{
    Intersections xs = (Intersections) malloc(sizeof(struct intersections));
    if (num > 0) {
        xs->xs = (Intersection) malloc(num * sizeof(struct intersection));
    } else {
        xs->xs = NULL;
    }

    xs->array_len = num;
    xs->num = 0;

    return xs;
}

void
intersection_free(Intersection x)
{
    if (x != NULL) {
        free(x);
    }
}

void
intersections_free(Intersections xs)
{
    if (xs != NULL) {
        if (xs->xs != NULL) {
            free(xs->xs);
        }
        free(xs);
    }
}

Ray
ray_array_alloc(double origin[4], double direction[4])
{
    Ray ray = (Ray) malloc(sizeof(struct ray));
    memcpy(ray->origin, origin, sizeof(ray->origin));
    memcpy(ray->direction, direction, sizeof(ray->direction));

    return ray;
}
Ray
ray_alloc(Point origin, Vector direction)
{
    return ray_array_alloc(origin->arr, direction->arr);
}

void
ray_free(Ray r) {
    if (r != NULL) {
        free(r);
    }
}

void
material(Material m)
{
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
}

Material
material_alloc()
{
    Material m = (Material) malloc(sizeof(struct material));
    material(m);
    return m;
}


void
ray_transform(Ray original, Matrix m, Ray res)
{
    matrix_array_multiply(m, original->origin, res->origin);
    matrix_array_multiply(m, original->direction, res->direction);
}

Ray
ray_transform_alloc(Ray original, Matrix m)
{
    Ray res = ray_array_alloc(original->origin, original->direction);
    ray_transform(original, m, res);
    return res;
}

int
ray_to_string(char *buf, size_t n, Ray r)
{
    return snprintf(buf, n, "Point: [%f %f %f] Vector: [%f %f %f]",
                    r->origin[0], r->origin[1], r->origin[2],
                    r->direction[0], r->direction[1], r->direction[2]);
}

Intersections
shape_intersect(Shape sh, Ray r)
{
    Matrix m = sh->transform_inverse;

    Ray r2 = ray_transform_alloc(r, m);

    Intersections xs = sh->local_intersect(sh, r2);

    ray_free(r2);

    return xs;
}

Vector
shape_normal_at(Shape sh, Point world_point, Intersection hit)
{
    Point local_point = sh->world_to_object(sh, world_point);
    Vector local_normal = sh->local_normal_at(sh, local_point, hit);
    Vector world_normal = sh->normal_to_world(sh, local_normal);
    Vector n = vector_normalize_alloc(world_normal);

    vector_free(world_normal);
    vector_free(local_normal);
    point_free(local_point);

    return n;
}

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

Vector
shape_normal_to_world(Shape sh, Vector local_normal)
{
    Matrix inv = sh->transform_inverse;
    Matrix tr = matrix_transpose_alloc(inv);
    
    Vector un_normal = matrix_vector_multiply_alloc(tr, local_normal);
    
    Vector normal = vector_normalize_alloc(un_normal);

    if (sh->parent != NULL) {
        Vector new_normal = sh->parent->normal_to_world(sh->parent, normal);
        vector_free(normal);
        normal = new_normal;
    }

    vector_free(un_normal);
    matrix_free(tr);

    return normal;
}

Point
shape_world_to_object(Shape sh, Point pt)
{
    if (sh->parent != NULL) {
        return sh->parent->world_to_object(sh->parent, pt);
    }

    Matrix m = sh->transform_inverse;
    Point p = matrix_point_multiply_alloc(m, pt);
    
    return p;
}

void
shape_divide(Shape sh, size_t threshold)
{
    return;
}

bool
shape_includes(Shape a, Shape b)
{
    return a == b;
}

void
shape_set_transform(Shape obj, Matrix m)
{
    if (obj != NULL) {
        if (obj->transform_inverse != NULL) {
            matrix_free(obj->transform_inverse);
        }
        if (obj->transform != NULL) {
            matrix_free(obj->transform);
        }
        obj->transform = m;
        obj->transform_inverse = matrix_inverse_alloc(m);
    }
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

int
shape_to_string(char *buf, size_t n, Shape sh)
{
    return snprintf(buf, n, "Shape:\n\ttransform: %p\n\tmaterial: %p\n\tparent: %p\n\ttype: %d\n",
        sh->transform,
        sh->material,
        sh->parent,
        sh->type);
}

