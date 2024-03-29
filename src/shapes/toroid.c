#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../libs/linalg/linalg.h"
#include "../libs/quartic/Roots3And4.h"

#include "shapes.h"
#include "bounding_box.h"

#include "toroid.h"


Intersections
toroid_local_intersect(Shape toroid, Ray r, bool stop_after_first_hit)
{
    double ox = r->origin[0];
    double oy = r->origin[1];
    double oz = r->origin[2];
    double dx = r->direction[0];
    double dy = r->direction[1];
    double dz = r->direction[2];

    double sum_d_sq = dx * dx + dy * dy + dz * dz;
    double e = ox * ox + oy * oy + oz * oz -
                toroid->fields.toroid.r1 * toroid->fields.toroid.r1 - 
                toroid->fields.toroid.r2 * toroid->fields.toroid.r2;
    double f = ox * dx + oy * dy + oz * dz;
    double four_a_sq = 4.0 * toroid->fields.toroid.r1 * toroid->fields.toroid.r1;
    double coeffs[5] = {
        e * e - four_a_sq * (toroid->fields.toroid.r2 * toroid->fields.toroid.r2 - oy * oy),
        4.0 * f * e + 2.0 * four_a_sq *oy * dy,
        2.0 * sum_d_sq * e + 4.0 * f * f + four_a_sq * dy * dy,
        4.0 * sum_d_sq * f,
        sum_d_sq * sum_d_sq
    };

    double solutions[4];
    int n = SolveQuartic(coeffs, solutions);
    if (n == 0) {
        return NULL;
    }

    Intersections xs = toroid->xs;
    Intersection x = xs->xs;
    xs->num = n;
    for (; n > 0; n--, x++) {
        intersection(solutions[n-1], toroid, x);
    }

    return xs;
}

void
toroid_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    double p_sq = sh->fields.toroid.r1 * sh->fields.toroid.r1 + sh->fields.toroid.r2 * sh->fields.toroid.r2;
    double mag = local_point[0] * local_point[0] + local_point[1] * local_point[1] + local_point[2] * local_point[2];
    Vector rv;
    rv[0] = 4.0 * local_point[0] * (mag - p_sq);
    rv[1] = 4.0 * local_point[1] * (mag - p_sq + 2.0 * sh->fields.toroid.r1 * sh->fields.toroid.r1);
    rv[2] = 4.0 * local_point[2] * (mag - p_sq);
    rv[3] = 0.0;
    vector_normalize(rv, res);
}

void
toroid_bounds(Shape toroid, Bounding_box *res)
{
    if (!toroid->bbox_valid) {
        toroid->bbox_valid = true;
        double arr[4] = {-toroid->fields.toroid.r1 - toroid->fields.toroid.r2, -toroid->fields.toroid.r2, -toroid->fields.toroid.r1 - toroid->fields.toroid.r2, 1.0};

        bounding_box_add_array(&(toroid->bbox), arr);

        arr[0] = toroid->fields.toroid.r1 + toroid->fields.toroid.r2;
        arr[1] = toroid->fields.toroid.r2;
        arr[2] = toroid->fields.toroid.r1 + toroid->fields.toroid.r2;
        bounding_box_add_array(&(toroid->bbox), arr);

        bounding_box_transform(&(toroid->bbox), toroid->transform, &(toroid->bbox_inverse));
    }

    *res = toroid->bbox;
}

void
toroid(Shape s)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_TOROID;
    s->fields.toroid.r1 = 0.75;
    s->fields.toroid.r2 = 0.25;

    bounding_box(&(s->bbox));
    s->bbox_valid = false;
    s->xs = intersections_empty(4);

    s->intersect = shape_intersect;
    s->local_intersect = toroid_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = toroid_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = toroid_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}

Shape
toroid_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    toroid(s);
    return s;
}
