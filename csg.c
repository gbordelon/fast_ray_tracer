#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "linalg.h"
#include "csg.h"
#include "shapes.h"
#include "renderer.h"

/*
 *
 * DIFFERENCE
 *
 * lhit | inr | inl | result |
 *    0 |   0 |   0 |      0 |
 *    0 |   0 |   1 |      1 |
 *    0 |   1 |   0 |      0 |
 *    0 |   1 |   1 |      1 |
 *    1 |   0 |   0 |      1 |
 *    1 |   0 |   1 |      1 |
 *    1 |   1 |   0 |      0 |
 *    1 |   1 |   1 |      0 |
 */

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

    for (i = 0, to = from = xs->xs; i < xs->num; i++, from++) {
        lhit = s->fields.csg.left->includes(s->fields.csg.left, from->object);

        if (intersection_allowed(s->fields.csg.op, lhit, inleft, inright)) {
            if (to != from) {
                *to = *from;
                //printf("\n");
            }
            num_remaining++;
            to++;
            //printf("keeping  x: %d %d %d\n", lhit, inleft, inright);
        } else {
            //printf("skipping x: %d %d %d\n", lhit, inleft, inright);
        }

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
    Bounding_box box = s->bounds(s);
    if (!bounding_box_intersects(box, r)) {
        return intersections_empty(0);
    }

    Intersections left_xs = s->fields.csg.left->intersect(s->fields.csg.left, r);
    Intersections right_xs = s->fields.csg.right->intersect(s->fields.csg.right, r);

    if ((left_xs->num + right_xs->num) > 0) {
        if (left_xs->num == 0) {
            // filter
            csg_filter_intersections(s, right_xs);
            intersections_free(left_xs);
            return right_xs;
        } else if (right_xs->num == 0) {
            // filter
            csg_filter_intersections(s, left_xs);
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

Bounding_box
csg_bounds(Shape csg)
{
    if (csg->bbox == NULL) {
        Bounding_box box = bounding_box_alloc();

        Bounding_box lbox = csg->fields.csg.left->parent_space_bounds(csg->fields.csg.left);
        Bounding_box rbox = csg->fields.csg.right->parent_space_bounds(csg->fields.csg.right);

        bounding_box_add_box(box, lbox);
        bounding_box_add_box(box, rbox);
        csg->bbox = box;
        csg->bbox_inverse =  bounding_box_transform(box, csg->transform);
    }
    return csg->bbox;
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
    s->bbox = NULL;
    s->bbox_inverse = NULL;

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

    s->bounds = csg_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}

Shape
csg_alloc(enum csg_ops_enum op, Shape left_child, Shape right_child)
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    csg(s, op, left_child, right_child);
    return s;
}
