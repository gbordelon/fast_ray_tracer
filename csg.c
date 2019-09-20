#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "linalg.h"
#include "csg.h"
#include "shapes.h"
#include "renderer.h"

bool
intersection_allowed(enum csg_ops_enum op, bool lhit, bool inl, bool inr)
{
    if (CSG_UNION == op) {
        return (lhit && (!inr)) || ((!lhit) && (!inl));
    } else if (CSG_INTERSECT == op) {
        return (lhit && inr) || ((!lhit) && inl);
    } else if (CSG_DIFFERENCE == op) {
        return (lhit && (!inr)) || ((!lhit) && inl);
    }

    printf("Unknown CSG operation.");
    return false;
}

void
csg_filter_intersections(Shape s, Intersections xs)
{
    int i;
    Intersection to, from;
    bool lhit;
    size_t num_remaining = 0;

    bool inleft = false;
    bool inright = false;

    for (i = 0, to = xs->xs, from = xs->xs; i < xs->num; i++) {
        lhit = s->fields.csg.left->includes(s->fields.csg.left, from->object);

        if (intersection_allowed(s->fields.csg.op, lhit, inleft, inright)) {
            if (to != from) {
                *to = *from;
            }
            num_remaining++;
            to++;
        }
        from++;

        if (lhit) {
            inleft = !inleft;
        } else {
            inright = !inright;
        }
    }
    xs->num = num_remaining;
}

Intersections
csg_local_intersect(Shape s, Ray r)
{
    Intersections left_xs = s->fields.csg.left->intersect(s->fields.csg.left, r);
    Intersections right_xs = s->fields.csg.right->intersect(s->fields.csg.right, r);

    if (left_xs->num + right_xs->num > 0) {
        if (left_xs->num == 0) {
            intersections_free(left_xs);
            return right_xs;
        } else if (right_xs->num == 0) {
            intersections_free(right_xs);
            return left_xs;
        }

        // both arrays together
        Intersections both = intersections_empty(left_xs->num + right_xs->num);
        memcpy(both->xs, left_xs->xs, left_xs->num * sizeof(struct intersection));
        memcpy(both->xs + left_xs->num, right_xs->xs, right_xs->num * sizeof(struct intersection));
        both->num = left_xs->num + right_xs->num;

        // sort both arrays together
        intersections_sort(both);

        // filter
        csg_filter_intersections(s, both);

        intersections_free(left_xs);
        intersections_free(right_xs);

        return both;
    }

    intersections_free(right_xs);

    return left_xs;
}

Vector
csg_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    printf("CSG local_normal_at called. This makes no sense. Returning null\n");
    return NULL;
}

bool
csg_includes(Shape csg, Shape other)
{
    return csg->fields.csg.left->includes(csg->fields.csg.left, other) ||
           csg->fields.csg.right->includes(csg->fields.csg.right, other);
}

void
csg_divide(Shape csg, size_t threshold)
{
    csg->fields.csg.left->divide(csg->fields.csg.left, threshold);
    csg->fields.csg.right->divide(csg->fields.csg.right, threshold);
}

void
csg(Shape s, enum csg_ops_enum op, Shape left_child, Shape right_child)
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_CSG;

    s->fields.csg.op = op;
    s->fields.csg.left = left_child;
    s->fields.csg.right = right_child;

    left_child->parent = s;
    right_child->parent = s;

    s->intersect = shape_intersect;
    s->local_intersect = csg_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = csg_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = csg_divide;
    s->includes = csg_includes;
}

Shape
csg_alloc(enum csg_ops_enum op, Shape left_child, Shape right_child)
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    csg(s, op, left_child, right_child);
    return s;
}
