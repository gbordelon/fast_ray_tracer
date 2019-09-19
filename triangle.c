#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "linalg.h"
#include "triangle.h"
#include "shapes.h"

Intersections
triangle_local_intersect(Shape triangle, Ray r)
{
    Vector dir_cross_e2 = vector_cross_arrays_alloc(r->direction, triangle->fields.triangle.e2);
    double det = array_dot(triangle->fields.triangle.e1, dir_cross_e2->arr);
    if (fabs(det) < EPSILON) {
        return intersections_empty(0);
    }

    double f = 1.0 / det;
    Vector p1_to_origin = vector_from_arrays_alloc(r->origin, triangle->fields.triangle.p1);
    double u = f * vector_dot(p1_to_origin, dir_cross_e2);
    if (u < 0 || u > 1) {
        return intersections_empty(0);
    }

    Vector origin_cross_e1 = vector_cross_arrays_alloc(p1_to_origin->arr, triangle->fields.triangle.e1);

    double v = f * array_dot(r->direction, origin_cross_e1->arr);
    if (v < 0 || ((u + v) > 1)) {
        return intersections_empty(0);
    }

    double t = f * array_dot(triangle->fields.triangle.e2, origin_cross_e1->arr);
    Intersections xs = intersections_empty(1);
    intersection(t, triangle, xs->xs);
    xs->num = 1;

    return xs;
}

Vector
triangle_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    return vector(sh->fields.triangle.u_normals.normal[0],
                  sh->fields.triangle.u_normals.normal[1], 
                  sh->fields.triangle.u_normals.normal[2]);
}


void
triangle(Shape s, double p1[4], double p2[4], double p3[4])
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_TRIANGLE;

    memcpy(s->fields.triangle.p1, p1, 4 * sizeof(double));
    memcpy(s->fields.triangle.p2, p2, 4 * sizeof(double));
    memcpy(s->fields.triangle.p3, p3, 4 * sizeof(double));

    array_from_arrays(p2, p1, &s->fields.triangle.e1[0]);
    array_from_arrays(p3, p1, &s->fields.triangle.e2[0]);

    Vector cross = vector_cross_arrays_alloc(s->fields.triangle.e2, s->fields.triangle.e1);
    Vector normal = vector_normalize_alloc(cross);
    memcpy(s->fields.triangle.u_normals.normal, normal->arr, sizeof(normal->arr));

    s->intersect = shape_intersect;
    s->local_intersect = triangle_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = triangle_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    vector_free(normal);
    vector_free(cross);
}


Shape
triangle_array_alloc(double p1[4], double p2[4], double p3[4])
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    triangle(s, p1, p2, p3);
    return s;
}

Shape
triangle_point_alloc(Point p1, Point p2, Point p3)
{
    return triangle_array_alloc(p1->arr, p2->arr, p3->arr);
}





Intersections
smooth_triangle_local_intersect(Shape triangle, Ray r)
{
    Vector dir_cross_e2 = vector_cross_arrays_alloc(r->direction, triangle->fields.triangle.e2);
    double det = array_dot(triangle->fields.triangle.e1, dir_cross_e2->arr);
    if (fabs(det) < EPSILON) {
        return intersections_empty(0);
    }

    double f = 1.0 / det;
    Vector p1_to_origin = vector_from_arrays_alloc(r->origin, triangle->fields.triangle.p1);
    double u = f * vector_dot(p1_to_origin, dir_cross_e2);
    if (u < 0 || u > 1) {
        return intersections_empty(0);
    }

    Vector origin_cross_e1 = vector_cross_arrays_alloc(p1_to_origin->arr, triangle->fields.triangle.e1);

    double v = f * array_dot(r->direction, origin_cross_e1->arr);
    if (v < 0 || ((u + v) > 1)) {
        return intersections_empty(0);
    }

    double t = f * array_dot(triangle->fields.triangle.e2, origin_cross_e1->arr);
    Intersections xs = intersections_empty(1);
    intersection_with_uv(t, u, v, triangle, xs->xs);
    xs->num = 1;

    return xs;
}

Vector
smooth_triangle_local_normal_at(Shape s, Point local_point, Intersection hit)
{
//         return self.n2 * hit.u + self.n3 * hit.v + self.n1 * (1 - hit.u - hit.v)
    Vector v2 = vector(s->fields.triangle.u_normals.s_normals.n2[0],
                       s->fields.triangle.u_normals.s_normals.n2[1],
                       s->fields.triangle.u_normals.s_normals.n2[2]);
    Vector v3 = vector(s->fields.triangle.u_normals.s_normals.n3[0],
                       s->fields.triangle.u_normals.s_normals.n3[1],
                       s->fields.triangle.u_normals.s_normals.n3[2]);
    Vector v1 = vector(s->fields.triangle.u_normals.s_normals.n1[0],
                       s->fields.triangle.u_normals.s_normals.n1[1],
                       s->fields.triangle.u_normals.s_normals.n1[2]);

    vector_scale(v2, hit->u);
    vector_scale(v3, hit->v);
    vector_scale(v1, 1 - hit->u - hit->v);

    v1->arr[0] += v2->arr[0] + v3->arr[0];
    v1->arr[1] += v2->arr[1] + v3->arr[1];
    v1->arr[2] += v2->arr[2] + v3->arr[2];

    vector_free(v2);
    vector_free(v3);

    return v1;
}


void
smooth_triangle(Shape s,
                double p1[4], double p2[4], double p3[4],
                double n1[4], double n2[4], double n3[4])
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_SMOOTH_TRIANGLE;

    memcpy(s->fields.triangle.p1, p1, 4 * sizeof(double));
    memcpy(s->fields.triangle.p2, p2, 4 * sizeof(double));
    memcpy(s->fields.triangle.p3, p3, 4 * sizeof(double));

    array_from_arrays(p2, p1, &s->fields.triangle.e1[0]);
    array_from_arrays(p3, p1, &s->fields.triangle.e2[0]);

    memcpy(s->fields.triangle.u_normals.s_normals.n1, n1, 4 * sizeof(double));
    memcpy(s->fields.triangle.u_normals.s_normals.n2, n2, 4 * sizeof(double));
    memcpy(s->fields.triangle.u_normals.s_normals.n3, n3, 4 * sizeof(double));

    s->intersect = shape_intersect;
    s->local_intersect = smooth_triangle_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = smooth_triangle_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;
}


Shape
smooth_triangle_array_alloc(double p1[4], double p2[4], double p3[4], double n1[4], double n2[4], double n3[4])
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    smooth_triangle(s, p1, p2, p3, n1, n2, n3);
    return s;
}

Shape
smooth_triangle_point_alloc(Point p1, Point p2, Point p3, Vector n1, Vector n2, Vector n3)
{
    return smooth_triangle_array_alloc(p1->arr, p2->arr, p3->arr, n1->arr, n2->arr, n3->arr);
}

