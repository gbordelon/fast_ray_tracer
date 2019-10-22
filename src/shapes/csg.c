#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "../libs/linalg/linalg.h"
#include "../renderer/renderer.h"

#include "csg.h"
#include "shapes.h"

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
            }
            num_remaining++;
            to++;
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
csg_local_intersect(Shape s, Ray r, bool stop_after_first_hit)
{
    Bounding_box box;
    s->bounds(s, &box);
    if (!bounding_box_intersects(&box, r)) {
        return NULL;
    }

    Intersections left_xs = s->fields.csg.left->intersect(s->fields.csg.left, r, stop_after_first_hit);
    Intersections right_xs = s->fields.csg.right->intersect(s->fields.csg.right, r, stop_after_first_hit);

    size_t total_num_intersections = 0;

    if (left_xs != NULL) {
        total_num_intersections += left_xs->num;
    }

    if (right_xs != NULL) {
        total_num_intersections += right_xs->num;
    }

    if (total_num_intersections > 0) {
        if (left_xs == NULL || left_xs->num == 0) {
            csg_filter_intersections(s, right_xs);
            return right_xs;
        } else if (right_xs == NULL || right_xs->num == 0) {
            csg_filter_intersections(s, left_xs);
            return left_xs;
        }

        if (s->xs->array_len <= left_xs->num + right_xs->num) {
            intersections_realloc(s->xs, left_xs->num + right_xs->num);
        }

        // both arrays together
        Intersections both = s->xs;
        both->num = 0;
        memcpy(both->xs, left_xs->xs, left_xs->num * sizeof(struct intersection));
        memcpy(both->xs + left_xs->num, right_xs->xs, right_xs->num * sizeof(struct intersection));
        both->num = left_xs->num + right_xs->num;

        // sort both arrays together
        intersections_sort(both);

        // filter
        csg_filter_intersections(s, both);

        return both;
    }

    return left_xs;
}

void
csg_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    printf("CSG local_normal_at called. This makes no sense.\n");
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
csg_bounds(Shape csg, Bounding_box *res)
{
    if (!csg->bbox_valid) {
        csg->bbox_valid = true;

        Bounding_box lbox, rbox;
        csg->fields.csg.left->parent_space_bounds(csg->fields.csg.left, &lbox);
        csg->fields.csg.right->parent_space_bounds(csg->fields.csg.right, &rbox);

        bounding_box_add_box(&(csg->bbox), &lbox);
        bounding_box_add_box(&(csg->bbox), &rbox);

        bounding_box_transform(&(csg->bbox), csg->transform, &(csg->bbox_inverse));
    }

    *res = csg->bbox;
}

void
csg(Shape s, enum csg_ops_enum op, Shape left_child, Shape right_child)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_CSG;
    bounding_box(&(s->bbox));
    s->bbox_valid = false;
    s->xs = intersections_empty(64);

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
