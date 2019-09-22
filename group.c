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
    Bounding_box box = group->bounds(group);
    if (!bounding_box_intersects(box, r)) {
        bounding_box_free(box);
        return intersections_empty(0);
    }

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
    bounding_box_free(box);

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
shape_swap(Shape l, Shape r)
{
    struct shape tmp = *l;
    *l = *r;
    *r = tmp;
}

struct partition_children_return_value {
    size_t left_count, middle_count, right_count;
    int left_start, middle_start, right_start;
};

struct partition_children_return_value
partition_children(Shape group)
{
    Bounding_box box = group->bounds(group);
    Bounding_box boxes = bounding_box_split_bounds(box);
    Bounding_box left = boxes;
    Bounding_box right = boxes+1;
    Bounding_box child_box;

    bool *left_map = (bool*)malloc(group->fields.group.num_children * sizeof(bool));
    bool *right_map = (bool*)malloc(group->fields.group.num_children * sizeof(bool));
    size_t left_count = 0, middle_count = 0, right_count = 0;
    int left_start = -1, middle_start = -1, right_start = -1;
    
    Shape from;
    int i, j;
    for (i = 0, from = group->fields.group.children;
            i < group->fields.group.num_children;
            i++, from++) {

        child_box = from->parent_space_bounds(from); // or just bounds?

        left_map[i] = false;
        right_map[i] = false;
        if (bounding_box_contains_box(left, child_box)) {
            left_map[i] = true;
            left_count++;
        } else if (bounding_box_contains_box(right, child_box)) {
            right_map[i] = true;
            right_count++;
        } else {
            middle_count++;
        }

        bounding_box_free(child_box);
    }

    // put left children into the beginning of the children array
    // then middle
    // then right
    // need from and to pointers initialized to 0 offset
    // while to and from are still pointing in the array:
    // if to is left, increment to and from
    // if to is not left, increment from until a left is found, then swap, then increment to and from
    // start another loop, to starts at left_count offset
    // while to and from are still pointing in the array:
    // if to is middle, increment to and from
    // if to is right, increment from until a middle is found, then swap, then increment to and from
    // array is now partitioned

    for (i = j = 0; i < group->fields.group.num_children && j < group->fields.group.num_children;) {
        if (left_map[i]) {
            if (left_start < 0) {
                left_start = i;
            }
            i++;
            j++;
        } else {
            while (j < group->fields.group.num_children && !left_map[j]) {
                j++;
            }
            if (j < group->fields.group.num_children) {
                shape_swap(group->fields.group.children + i, group->fields.group.children + j);
            }
        }
    }

    // i should equal left_count
    // if left_count is zero then left_start should equal -1
    // else left_start should equal 0

    for (j = i; i < group->fields.group.num_children && j < group->fields.group.num_children;) {
        if (!right_map[i]) { // so middle
            if (middle_start < 0) {
                middle_start = i;
            }
            i++;
            j++;
        } else {
            while (j < group->fields.group.num_children && right_map[j]) {
                j++;
            }
            if (j < group->fields.group.num_children) {
                shape_swap(group->fields.group.children + i, group->fields.group.children + j);
            }
        }
    }

    // i should equal left_count + middle_count
    // if middle_count is zero then middle_start should equal -1
    // else middle_start should equal left_count
    // left_count + middle_count + right_count should equal group->fields.group.num_children

    right_start = i;

    struct partition_children_return_value rv;
    rv.left_count = left_count;
    rv.middle_count = middle_count;
    rv.right_count = right_count;
    rv.left_start = left_start;
    rv.middle_start = middle_start;
    rv.right_start = right_start;

    free(right_map);
    free(left_map);
    bounding_box_free(boxes);
    bounding_box_free(box);

    return rv;
}

void
group_divide(Shape g, size_t threshold)
{
    if (threshold < g->fields.group.num_children) {
        // partition children
        struct partition_children_return_value map = partition_children(g);
        Shape new_children = NULL, new_group_pos;
        size_t new_array_size = map.middle_count;

        // make subgroup for each partition
        if (map.middle_count != g->fields.group.num_children) {
            if (map.left_count > 0) {
                new_array_size++;
            }
            if (map.right_count > 0) {
                new_array_size++;
            }
            new_children = (Shape) malloc(new_array_size * sizeof(struct shape));
            new_group_pos = new_children;
        }

        if (map.left_count > 0) {
            // make a sub group from left_start using left_count as num_children
            Shape left_children = (Shape) malloc(map.left_count * sizeof(struct shape));
            memcpy(left_children, g->fields.group.children + map.left_start, map.left_count * sizeof(struct shape));
            group(new_group_pos, left_children, map.left_count);
            new_group_pos->fields.group.children_need_free = true;
            new_group_pos->parent = g;
            new_group_pos++;
        }
        if (map.right_count > 0) {
            // make a sub group from right_start using right_count as num_children
            Shape right_children = (Shape)malloc(map.right_count * sizeof(struct shape));
            memcpy(right_children, g->fields.group.children + map.right_start, map.right_count * sizeof(struct shape));
            group(new_group_pos, right_children, map.right_count);
            new_group_pos->fields.group.children_need_free = true;
            new_group_pos->parent = g;
            new_group_pos++;
        }
        if (map.middle_count != g->fields.group.num_children) {
            // copy from middle_start for middle_count into new array
            memcpy(new_group_pos, g->fields.group.children + map.middle_start, map.middle_count * sizeof(struct shape));
            if (g->fields.group.children_need_free) {
                free(g->fields.group.children);
            }
            g->fields.group.children = new_children;
            g->fields.group.children_need_free = true;
            g->fields.group.num_children = new_array_size;
        }
    }
    
    Shape child;
    int i;
    for (i = 0, child = g->fields.group.children;
            i < g->fields.group.num_children;
            i++, child++) {
        child->divide(child, threshold);
    }
}

Bounding_box
group_bounds_alloc(Shape group)
{
    Bounding_box box = bounding_box_alloc();
    Shape child;
    Bounding_box bbox;
    int i;
    for (i = 0, child = group->fields.group.children;
            i < group->fields.group.num_children;
            i++, child++) {
        bbox = child->parent_space_bounds(child);
        bounding_box_add_box(box, bbox);
        bounding_box_free(bbox);
    }

    return box;
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

    s->bounds = group_bounds_alloc;
    s->parent_space_bounds = shape_parent_space_bounds_alloc;
}

Shape
group_alloc(Shape children, size_t num_children)
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    group(s, children, num_children);
    return s;
}
