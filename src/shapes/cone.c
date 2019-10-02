#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "../libs/linalg/linalg.h"

#include "cone.h"
#include "shapes.h"

bool
check_cap(Ray r, double t, double y)
{
    double x = r->origin[0] + t * r->direction[0];
    double z = r->origin[2] + t * r->direction[2];
    double sm = x * x + z * z;
    return sm <= fabs(y);
}

void
intersect_caps(Shape cone, Ray r, Intersections xs)
{
    if ((!cone->fields.cone.closed) || equal(r->direction[1], 0.0)) {
        return;
    }

    Intersection x = xs->xs + xs->num;

    double t = (cone->fields.cone.minimum - r->origin[1]) / r->direction[1];
    if (check_cap(r, t, cone->fields.cone.minimum)) {
        intersection(t, cone, x++);
        xs->num++;
    }

    t = (cone->fields.cone.maximum - r->origin[1]) / r->direction[1];
    if (check_cap(r, t, cone->fields.cone.maximum)) {
        intersection(t, cone, x);
        xs->num++;
    }
}

Intersections
cone_local_intersect(Shape cone, Ray r)
{
    Intersections xs = cone->xs;
    xs->num = 0;

    Intersection x = xs->xs;

    double a = r->direction[0] * r->direction[0] +
               r->direction[2] * r->direction[2] -
               r->direction[1] * r->direction[1];
    double b = 2 * (r->origin[0] * r->direction[0] +
               r->origin[2] * r->direction[2] -
               r->origin[1] * r->direction[1]);
    double c = r->origin[0] * r->origin[0] +
               r->origin[2] * r->origin[2] -
               r->origin[1] * r->origin[1];

    if (equal(a,0.0)) {
        if (!equal(b,0.0)) {
            intersection(-c / (2 * b), cone, x++);
            xs->num++;
        }
    } else {
        double disc = b * b - 4 * a * c;
        if (disc < 0) {
            return xs;
        }

        double discsqrt = sqrt(disc);
        double t0 = (-b - discsqrt) / (2 * a);
        double t1 = (-b + discsqrt) / (2 * a);
        if (t0 > t1) {
            double tmp = t0;
            t0 = t1;
            t1 = tmp;
        }

        double y0 = r->origin[1] + t0 * r->direction[1];
        if (cone->fields.cone.minimum < y0 && y0 < cone->fields.cone.maximum) {
            intersection(t0, cone, x++);
            xs->num++;
        }

        double y1 = r->origin[1] + t1 * r->direction[1];
        if (cone->fields.cone.minimum < y1 && y1 < cone->fields.cone.maximum) {
            intersection(t1, cone, x++);
            xs->num++;
        }
    }

    intersect_caps(cone, r, xs);

    return xs;
}

void
cone_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    double dist = local_point[0] * local_point[0] + local_point[2] * local_point[2];
    vector_default(res);

    if (dist < 1 && ((sh->fields.cone.maximum - EPSILON) <= local_point[1])) {
        res[1] = 1;
    } else if (dist < 1 && ((sh->fields.cone.minimum + EPSILON) >= local_point[1])) {
        res[1] = -1;
    } else {
        double y = sqrt(dist);
        if (local_point[1] > 0) {
            y = -y;
        }

        res[0] = local_point[0];
        res[1] = y;
        res[2] = local_point[2];
    }
}

Bounding_box
cone_bounds(Shape cone)
{
    if (cone->bbox == NULL) {
        double a = fabs(cone->fields.cone.minimum);
        double b = fabs(cone->fields.cone.maximum);
        double limit = fmax(a,b);
        double arr[4] = {-limit, cone->fields.cone.minimum, -limit, 1.0};
        Bounding_box box = bounding_box_alloc();

        bounding_box_add_array(box, arr);

        arr[0] = limit;
        arr[1] = cone->fields.cone.maximum;
        arr[2] = limit;
        bounding_box_add_array(box, arr);

        cone->bbox = box;
        cone->bbox_inverse = bounding_box_transform(box, cone->transform);
    }

    return cone->bbox;
}

void
cone(Shape s)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_CONE;
    s->bbox = NULL;
    s->bbox_inverse = NULL;
    s->xs = intersections_empty(4);

    s->fields.cone.minimum = -DBL_MAX;
    s->fields.cone.maximum = DBL_MAX;
    s->fields.cone.closed = false;

    s->intersect = shape_intersect;
    s->local_intersect = cone_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = cone_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = cone_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}

Shape
cone_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    cone(s);
    return s;
}
