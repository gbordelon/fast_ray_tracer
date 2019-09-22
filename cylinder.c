#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "linalg.h"
#include "cylinder.h"
#include "shapes.h"

bool
cylinder_check_cap(Ray r, double t)
{
    double x = r->origin[0] + t * r->direction[0];
    double z = r->origin[2] + t * r->direction[2];
    double sm = x * x + z * z;
    return sm <= 1;
}

void
cylinder_intersect_caps(Shape cylinder, Ray r, Intersections xs)
{
    if ((!cylinder->fields.cylinder.closed) || equal(r->direction[1], 0.0)) {
        return;
    }

    Intersection x = xs->xs + xs->num;
    double t1 = (cylinder->fields.cylinder.minimum - r->origin[1]) / r->direction[1];
    double t2 = (cylinder->fields.cylinder.maximum - r->origin[1]) / r->direction[1];

    if (cylinder_check_cap(r, t1)) {
        intersection(t1, cylinder, x++);
        xs->num++;
    }

    if (cylinder_check_cap(r, t2)) {
        intersection(t2, cylinder, x);
        xs->num++;
    }
}

Intersections
cylinder_local_intersect(Shape cylinder, Ray r)
{
    Intersections xs = intersections_empty(4);
    Intersection x = xs->xs;

    double a = r->direction[0] * r->direction[0] +
               r->direction[2] * r->direction[2];
    double b = 2 * (r->origin[0] * r->direction[0] +
               r->origin[2] * r->direction[2]);
    double c = r->origin[0] * r->origin[0] +
               r->origin[2] * r->origin[2] -
               1;

    if (!equal(a,0.0)) {
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
        if (cylinder->fields.cylinder.minimum <= y0 && y0 <= cylinder->fields.cylinder.maximum) {
            intersection(t0, cylinder, x++);
            xs->num++;
        }

        double y1 = r->origin[1] + t1 * r->direction[1];
        if (cylinder->fields.cylinder.minimum <= y1 && y1 <= cylinder->fields.cylinder.maximum) {
            intersection(t1, cylinder, x);
            xs->num++;
        }
    }

    cylinder_intersect_caps(cylinder, r, xs);

    return xs;
}

Vector
cylinder_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    double dist = local_point->arr[0] * local_point->arr[0] +
                  local_point->arr[2] * local_point->arr[2];
    if (dist < 1 && ((sh->fields.cylinder.maximum - EPSILON) <= local_point->arr[1])) {
        return vector(0,1,0);
    } else if (dist < 1 && ((sh->fields.cylinder.minimum + EPSILON) >= local_point->arr[1])) {
        return vector(0,-1,0);
    }

    return vector(local_point->arr[0], 0, local_point->arr[2]);
}

Bounding_box
cylinder_bounds_alloc(Shape cylinder)
{
    double arr[4] = {-1.0, cylinder->fields.cylinder.minimum, -1.0, 1.0};
    Bounding_box box = bounding_box_alloc();

    bounding_box_add_array(box, arr);

    arr[0] = 1;
    arr[1] = cylinder->fields.cylinder.maximum;
    arr[2] = 1;
    bounding_box_add_array(box, arr);

    return box;
}


void
cylinder(Shape s)
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_CYLINDER;

    s->fields.cylinder.minimum = -INFINITY;
    s->fields.cylinder.maximum = INFINITY;
    s->fields.cylinder.closed = false;

    s->intersect = shape_intersect;
    s->local_intersect = cylinder_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = cylinder_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = shape_divide;
    s->includes = shape_includes;

    s->bounds = cylinder_bounds_alloc;
    s->parent_space_bounds = shape_parent_space_bounds_alloc;
}

Shape
cylinder_alloc()
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    cylinder(s);
    return s;
}
