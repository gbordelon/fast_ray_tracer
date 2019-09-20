#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "linalg.h"
#include "group.h"
#include "shapes.h"
#include "renderer.h"



void
group_add_children(Shape group, Shape children, size_t num_children)
{
    if (group == children) {
        printf("Error trying to add a group to its own children list\n");
        return;
    }

    // filter out dupes so only only of each child appears in the child array
    Shape child, new_child;
    int i, j;
    size_t num_children_to_add = num_children;
    unsigned char *truth_table = (unsigned char*) malloc(num_children * sizeof(unsigned char));
    for (i = 0, child = group->fields.group.children;
            i < group->fields.group.num_children;
            i++, child++) {
        for (j = 0, new_child = children;
                j < num_children;
                j++, new_child++) {
            *(truth_table + j) = 1;
            if (new_child == child) {
                *(truth_table + j) = 0;
                num_children_to_add--;
            }
        }
    }

    if (num_children_to_add > 0) {
        Shape total_children = (Shape)malloc((num_children_to_add + group->fields.group.num_children) * sizeof(struct shape));
        memcpy(total_children, group->fields.group.children, group->fields.group.num_children * sizeof(struct shape));

        for (j = 0, new_child = children, child = total_children + group->fields.group.num_children;
                j < num_children;
                j++, new_child++) {
            if (*(truth_table + j) == 1) {
                *child = *new_child;
                child->parent = group;
                child++;
            }
        }
        if (group->fields.group.children_need_free) {
            free(group->fields.group.children);
        }
        group->fields.group.children_need_free = true;
        group->fields.group.children = total_children;
        group->fields.group.num_children += num_children_to_add;
    }

    free(truth_table);
}

Intersections
group_local_intersect(Shape group, Ray r)
{
    Intersections *children_xs = (Intersections *) malloc(group->fields.group.num_children * sizeof(Intersection));

    int i;
    Shape child;
    size_t total_xs_num = 0;
    for (i = 0, child = group->fields.group.children;
            i < group->fields.group.num_children;
            i++, child++) {
        *(children_xs + i) = child->intersect(child, r);
        total_xs_num += (*(children_xs + i))->num;
    }

    if (total_xs_num == 0) {
        for (i = 0; i < group->fields.group.num_children; i++) {
            intersections_free(*(children_xs + i));
        }

        free(children_xs);
        return intersections_empty(0);
    }

    Intersections all_xs = intersections_empty(total_xs_num);

    for (i = 0; i < group->fields.group.num_children; i++) {
        memcpy(all_xs->xs + all_xs->num,
               (*(children_xs + i))->xs,
               (*(children_xs + i))->num * sizeof(struct intersection));
        all_xs->num += (*(children_xs + i))->num;
        intersections_free(*(children_xs + i));
    }

    intersections_sort(all_xs);

    free(children_xs);

    return all_xs;
}

Vector
group_local_normal_at(Shape sh, Point local_point, Intersection hit)
{
    printf("Group local_normal_at called. This makes no sense. Returning null\n");
    return NULL;
}

bool
group_includes(Shape group, Shape other)
{
    Shape child;
    int i;
    for (i = 0, child = group->fields.group.children;
            i < group->fields.group.num_children;
            i++, child++) {
        if (child->includes(child, other)) {
            return true;
        }
    }
    return false;
}

void
group_divide(Shape group, size_t threshold)
{
    if (threshold < group->fields.group.num_children) {
        // partition children
        // make subgroup for each partition
    }
    
    Shape child;
    int i;
    for (i = 0, child = group->fields.group.children;
            i < group->fields.group.num_children;
            i++, child++) {
        child->divide(child, threshold);
    }
}

void
group(Shape s, Shape children, size_t num_children)
{
    s->transform = NULL;
    s->transform_inverse = NULL;

    shape_set_transform(s, matrix_identity_alloc());

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_GROUP;

    s->fields.group.children = children;
    s->fields.group.num_children = num_children;
    s->fields.group.children_need_free = false;

    Shape child;
    int i;
    for (i = 0, child = children; i < num_children; i++, child++) {
        child->parent = s;
    }

    s->intersect = shape_intersect;
    s->local_intersect = group_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = group_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = group_divide;
    s->includes = group_includes;
}

Shape
group_alloc(Shape children, size_t num_children)
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    group(s, children, num_children);
    return s;
}
