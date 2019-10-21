#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../libs/linalg/linalg.h"

#include "triangle.h"
#include "shapes.h"

Intersections
triangle_local_intersect(Shape triangle, Ray r)
{
    Vector dir_cross_e2;
    vector_cross(r->direction, triangle->fields.triangle.e2, dir_cross_e2);
    double det = vector_dot(triangle->fields.triangle.e1, dir_cross_e2);
    if (fabs(det) < EPSILON) {
        return NULL;
    }

    double f = 1.0 / det;
    Vector p1_to_origin;
    vector_from_points(r->origin, triangle->fields.triangle.p1, p1_to_origin);
    double u = f * vector_dot(p1_to_origin, dir_cross_e2);
    if (u < 0 || u > 1) {
        return NULL;
    }

    Vector origin_cross_e1;
    vector_cross(p1_to_origin, triangle->fields.triangle.e1, origin_cross_e1);
    double v = f * vector_dot(r->direction, origin_cross_e1);
    if (v < 0 || ((u + v) > 1)) {
        return NULL;
    }

    double t = f * vector_dot(triangle->fields.triangle.e2, origin_cross_e1);
    Intersections xs = triangle->xs;
    xs->num = 0;

    //intersection(t, triangle, xs->xs);
    intersection_with_uv(t, u, v, triangle, xs->xs);
    xs->num = 1;

    return xs;
}

void
triangle_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    memcpy(res, sh->fields.triangle.u_normals.normal, sizeof(Vector));
}

void
triangle_bounds(Shape triangle, Bounding_box *res)
{
    if (!triangle->bbox_valid) {
        triangle->bbox_valid = true;

        bounding_box_add_array(&(triangle->bbox), triangle->fields.triangle.p1);
        bounding_box_add_array(&(triangle->bbox), triangle->fields.triangle.p2);
        bounding_box_add_array(&(triangle->bbox), triangle->fields.triangle.p3);
        bounding_box_transform(&(triangle->bbox), triangle->transform, &(triangle->bbox_inverse));
    }

    *res = triangle->bbox;
}

void
triangle(Shape s, Point p1, Point p2, Point p3)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_TRIANGLE;
    bounding_box(&(s->bbox));
    s->bbox_valid = false;
    s->xs = intersections_empty(1);
    s->fields.triangle.use_textures = false;

    memcpy(s->fields.triangle.p1, p1, sizeof(Point));
    memcpy(s->fields.triangle.p2, p2, sizeof(Point));
    memcpy(s->fields.triangle.p3, p3, sizeof(Point));

    vector_from_points(p2, p1, s->fields.triangle.e1);
    vector_from_points(p3, p1, s->fields.triangle.e2);

    Vector cross;
    vector_cross(s->fields.triangle.e2, s->fields.triangle.e1, cross);
    Vector normal;
    vector_normalize(cross, normal);
    memcpy(s->fields.triangle.u_normals.normal, normal, sizeof(Vector));

    s->intersect = shape_intersect;
    s->local_intersect = triangle_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = triangle_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = triangle_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}


Shape
triangle_array_alloc(Point p1, Point p2, Point p3)
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    triangle(s, p1, p2, p3);
    return s;
}

Shape
triangle_point_alloc(Point p1, Point p2, Point p3)
{
    return triangle_array_alloc(p1, p2, p3);
}

Intersections
smooth_triangle_local_intersect(Shape triangle, Ray r)
{
    Vector dir_cross_e2;
    vector_cross(r->direction, triangle->fields.triangle.e2, dir_cross_e2);
    double det = vector_dot(triangle->fields.triangle.e1, dir_cross_e2);
    if (fabs(det) < EPSILON) {
        return NULL;
    }

    double f = 1.0 / det;
    Vector p1_to_origin;
    vector_from_points(r->origin, triangle->fields.triangle.p1, p1_to_origin);
    double u = f * vector_dot(p1_to_origin, dir_cross_e2);
    if (u < 0 || u > 1) {
        return NULL;
    }

    Vector origin_cross_e1;
    vector_cross(p1_to_origin, triangle->fields.triangle.e1, origin_cross_e1);

    double v = f * vector_dot(r->direction, origin_cross_e1);
    if (v < 0 || ((u + v) > 1)) {
        return NULL;
    }

    double t = f * vector_dot(triangle->fields.triangle.e2, origin_cross_e1);
    Intersections xs = triangle->xs;
    xs->num = 0;

    intersection_with_uv(t, u, v, triangle, xs->xs);
    xs->num = 1;

    return xs;
}

void
smooth_triangle_local_normal_at(Shape s, Point local_point, Intersection hit, Vector res)
{
    Vector v2;
    Vector v3;
    vector_copy(res, s->fields.triangle.u_normals.s_normals.n1);
    vector_copy(v2, s->fields.triangle.u_normals.s_normals.n2);
    vector_copy(v3, s->fields.triangle.u_normals.s_normals.n3);

    vector_scale(v2, hit->u);
    vector_scale(v3, hit->v);
    vector_scale(res, 1.0 - hit->u - hit->v);

    res[0] += v2[0] + v3[0];
    res[1] += v2[1] + v3[1];
    res[2] += v2[2] + v3[2];

}


void
smooth_triangle(Shape s,
                Point p1, Point p2, Point p3,
                Vector n1, Vector n2, Vector n3)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_SMOOTH_TRIANGLE;
    bounding_box(&(s->bbox));
    s->bbox_valid = false;
    s->xs = intersections_empty(1);
    s->fields.triangle.use_textures = false;

    memcpy(s->fields.triangle.p1, p1, sizeof(Point));
    memcpy(s->fields.triangle.p2, p2, sizeof(Point));
    memcpy(s->fields.triangle.p3, p3, sizeof(Point));

    vector_from_points(p2, p1, s->fields.triangle.e1);
    vector_from_points(p3, p1, s->fields.triangle.e2);

    vector_copy(s->fields.triangle.u_normals.s_normals.n1, n1);
    vector_copy(s->fields.triangle.u_normals.s_normals.n2, n2);
    vector_copy(s->fields.triangle.u_normals.s_normals.n3, n3);

    s->intersect = shape_intersect;
    s->local_intersect = smooth_triangle_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = smooth_triangle_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = triangle_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}

void
convert_triangle_to_smooth_triangle(Shape s, Vector n1, Vector n2, Vector n3)
{
    if (s->type != SHAPE_TRIANGLE) {
        return;
    }
    s->type = SHAPE_SMOOTH_TRIANGLE;
    vector_copy(s->fields.triangle.u_normals.s_normals.n1, n1);
    vector_copy(s->fields.triangle.u_normals.s_normals.n2, n3);
    vector_copy(s->fields.triangle.u_normals.s_normals.n3, n3);
    s->local_intersect = smooth_triangle_local_intersect;
    s->local_normal_at = smooth_triangle_local_normal_at;
}
