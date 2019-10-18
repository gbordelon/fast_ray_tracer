#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "../libs/linalg/linalg.h"
#include "../renderer/renderer.h"
#include "../renderer/world.h"

#include "shapes.h"
#include "group.h"

void
recursive_print(Shape sh, size_t indent)
{
    Shape child;
    int i;

    char *spaces = (char *)malloc((indent + 1) * sizeof(char));
    for (i = 0; i < indent; ++i) {
        *(spaces + i) = ' ';
    }
    *(spaces+i) = '\0';

    printf("%s", spaces);
    if (sh->type == SHAPE_GROUP) {
        printf("%lu\n", sh->fields.group.num_children);
    } else if (sh->type != SHAPE_TRIANGLE) {
        printf("type: %d\n", sh->type);
    }

    if (sh->type == SHAPE_GROUP) {
        for (i = 0, child = sh->fields.group.children;
                i < sh->fields.group.num_children;
                i++, child++) {
            if (child->type == SHAPE_GROUP) {
                recursive_print(child, indent + 2);
            }        }
    } else if (sh->type == SHAPE_CSG) {
        recursive_print(sh->fields.csg.left, indent + 2);
        recursive_print(sh->fields.csg.right, indent + 2);
    }

    free(spaces);
}


void
group_recursive_parent_update(Shape sh, Shape parent)
{
    Shape child;
    int i;

    sh->parent = parent;

    if (sh->type == SHAPE_GROUP) {
        for (i = 0, child = sh->fields.group.children;
                i < sh->fields.group.num_children;
                i++, child++) {
            group_recursive_parent_update(child, sh);
        }
    } else if (sh->type == SHAPE_CSG) {
        group_recursive_parent_update(sh->fields.csg.left, sh);
        group_recursive_parent_update(sh->fields.csg.right, sh);
    }
}

void
invalidate_bounding_box(Shape group)
{
    bounding_box_free(group->bbox);
    bounding_box_free(group->bbox_inverse);

    group->bbox = NULL;
    group->bbox_inverse = NULL;
}

void
recursive_invalidate_bounding_box(Shape sh)
{
    Shape child;
    int i;

    invalidate_bounding_box(sh);

    if (sh->type == SHAPE_GROUP) {
        for (i = 0, child = sh->fields.group.children;
                i < sh->fields.group.num_children;
                i++, child++) {
            recursive_invalidate_bounding_box(child);
        }
    } else if (sh->type == SHAPE_CSG) {
        recursive_invalidate_bounding_box(sh->fields.csg.left);
        recursive_invalidate_bounding_box(sh->fields.csg.right);
    }
}

void
group_add_children_stage(Shape group, Shape children, size_t num_children)
{
    int from;
    for (from = 0; from < num_children; from++) {
        // allow dupes but not self reference
        if ((children + from) != group) {
            // if there is not enough room, get more
            if (group->fields.group.num_children >= (group->fields.group.size_children_array - 1)) {
                group->fields.group.children = (Shape) realloc(group->fields.group.children, group->fields.group.size_children_array * 2 * sizeof(struct shape));

                group->fields.group.size_children_array *= 2;
            }

            // copy child_from to child_to
            shape_copy(children + from, group, group->fields.group.children + group->fields.group.num_children);

            // increment counts and pointers
            group->fields.group.num_children += 1;
        }
    }
}

void
group_add_children_finish(Shape group)
{
    if (group != NULL) {
        recursive_invalidate_bounding_box(group);
        group_recursive_parent_update(group, group->parent);
        if (group->fields.group.num_children > 0) {
            group->fields.group.children_xs = (Intersections *) realloc(group->fields.group.children_xs, group->fields.group.num_children * sizeof(Intersections));
        }
    }
}

void
group_add_children(Shape group, Shape children, size_t num_children)
{
    group_add_children_stage(group, children, num_children);
    group_add_children_finish(group);
}

Intersections
group_local_intersect(Shape group, Ray r)
{
    Bounding_box box = group->bounds(group);
    if (!bounding_box_intersects(box, r)) {
        return NULL;
    }

    Intersections *children_xs = group->fields.group.children_xs;

    int i;
    Shape child;
    size_t total_xs_num = 0;
    for (i = 0, child = group->fields.group.children;
            i < group->fields.group.num_children;
            i++, child++) {
        *(children_xs + i) = child->intersect(child, r);
        if (*(children_xs + i) != NULL) {
            total_xs_num += (*(children_xs + i))->num;
        }
    }

    if (total_xs_num == 0) {
        return NULL;
    }

    Intersections all_xs = group->xs;
    all_xs->num = 0;
    if (all_xs->array_len <= total_xs_num) {
        intersections_realloc(all_xs, total_xs_num);
    }

    for (i = 0; i < group->fields.group.num_children; i++) {
        if (*(children_xs + i) != NULL) {
            memcpy(all_xs->xs + all_xs->num,
                   (*(children_xs + i))->xs,
                   (*(children_xs + i))->num * sizeof(struct intersection));
            all_xs->num += (*(children_xs + i))->num;
        }
    }

    intersections_sort(all_xs);

    return all_xs;
}

void
group_local_normal_at(Shape sh, Point local_point, Intersection hit, Vector res)
{
    printf("Group local_normal_at called. This makes no sense. Returning null\n");
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

        child_box = from->parent_space_bounds(from);

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
    }

    // put left children into the beginning of the children array
    // then middle
    // then right
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
                *(left_map + i) ^= *(left_map + j);
                *(left_map + j) ^= *(left_map + i);
                *(left_map + i) ^= *(left_map + j);
                *(right_map + i) ^= *(right_map + j);
                *(right_map + j) ^= *(right_map + i);
                *(right_map + i) ^= *(right_map + j);
            }
        }
    }

    // i should equal left_count
    // if left_count is zero then left_start should equal -1
    // else left_start should equal 0
    for (j = i; i < group->fields.group.num_children && j < group->fields.group.num_children;) {
        if (!right_map[i]) {
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
                *(left_map + i) ^= *(left_map + j);
                *(left_map + j) ^= *(left_map + i);
                *(left_map + i) ^= *(left_map + j);
                *(right_map + i) ^= *(right_map + j);
                *(right_map + j) ^= *(right_map + i);
                *(right_map + i) ^= *(right_map + j);
            }
        }
    }

    // i should equal left_count + middle_count
    // if middle_count is zero then middle_start should equal -1
    // else middle_start should equal left_count
    // left_count + middle_count + right_count should equal group->fields.group.num_children
    if (i < group->fields.group.num_children) {
        right_start = i;
    }

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

    return rv;
}

void
group_divide(Shape g, size_t threshold)
{
    Shape child;
    int i;

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
            if (new_array_size < g->fields.group.size_children_array) {
                new_array_size = g->fields.group.size_children_array;
            }
            new_children = array_of_shapes(new_array_size);
            new_group_pos = new_children;
        }

        if (map.left_count > 0) {
            // make a sub group from left_start using left_count as num_children
            group(new_group_pos, g->fields.group.children + map.left_start, map.left_count);
            for (i = map.left_start; i < map.left_start + map.left_count; ++i) {
                shape_free(g->fields.group.children + i);
            }
            // recursive parent update on children
            group_recursive_parent_update(new_group_pos, g);
            new_group_pos++;
        }
        if (map.right_count > 0) {
            // make a sub group from right_start using right_count as num_children
            group(new_group_pos, g->fields.group.children + map.right_start, map.right_count);
            for (i = map.right_start; i < map.right_start + map.right_count; ++i) {
                shape_free(g->fields.group.children + i);
            }
            // recursive parent update on children
            group_recursive_parent_update(new_group_pos, g);
            new_group_pos++;
        }
        if (map.middle_count != g->fields.group.num_children) {
            // copy from middle_start for middle_count into new array
            memcpy(new_group_pos, g->fields.group.children + map.middle_start, map.middle_count * sizeof(struct shape));
            
            free(g->fields.group.children);
            g->fields.group.children = new_children;
            g->fields.group.num_children = map.middle_count;
            if (map.left_count > 0) {
                g->fields.group.num_children++;
            }
            if (map.right_count > 0) {
                g->fields.group.num_children++;
            }
            g->fields.group.children_xs = (Intersections *) realloc(g->fields.group.children_xs, g->fields.group.num_children * sizeof(Intersections));
            g->fields.group.size_children_array = new_array_size;
            group_recursive_parent_update(g, g->parent);
        }
    }
    
    for (i = 0, child = g->fields.group.children;
            i < g->fields.group.num_children;
            i++, child++) {
        child->divide(child, threshold);
    }
}

Bounding_box
group_bounds(Shape group)
{
    if (group->bbox == NULL) {
        Bounding_box box = bounding_box_alloc();
        Shape child;
        Bounding_box bbox;
        int i;
        for (i = 0, child = group->fields.group.children;
                i < group->fields.group.num_children;
                i++, child++) {
            bbox = child->parent_space_bounds(child);
            bounding_box_add_box(box, bbox);
        }
        group->bbox = box;
        group->bbox_inverse = bounding_box_transform(box, group->transform);
    }
    return group->bbox;
}

#define DEFAULT_MIN_CHILD_LEN 16
void
group(Shape s, Shape children, size_t num_children)
{
    shape_set_transform(s, MATRIX_IDENTITY);

    s->material = material_alloc();
    s->parent = NULL;
    s->type = SHAPE_GROUP;
    s->xs = intersections_empty(64);

    size_t array_len = DEFAULT_MIN_CHILD_LEN > num_children ? DEFAULT_MIN_CHILD_LEN : num_children;
    s->fields.group.children = (Shape) malloc(array_len * sizeof(struct shape));
    if (num_children > 0) {
        int i;
        for (i = 0; i < num_children; ++i) {
            shape_copy(children + i, s, s->fields.group.children + i);
        }
        s->fields.group.children_xs = (Intersections *) malloc(num_children * sizeof(Intersections));
    } else {
        s->fields.group.children_xs = NULL;
    }

    s->fields.group.num_children = num_children;
    s->fields.group.size_children_array = array_len;
    s->bbox = NULL;
    s->bbox_inverse = NULL;

    group_recursive_parent_update(s, s->parent);
    recursive_invalidate_bounding_box(s);

    s->intersect = shape_intersect;
    s->local_intersect = group_local_intersect;
    s->normal_at = shape_normal_at;
    s->local_normal_at = group_local_normal_at;
    s->normal_to_world = shape_normal_to_world;
    s->world_to_object = shape_world_to_object;
    s->divide = group_divide;
    s->includes = group_includes;

    s->bounds = group_bounds;
    s->parent_space_bounds = shape_parent_space_bounds;
}

Shape
group_alloc(Shape children, size_t num_children)
{
    Shape s = (Shape) malloc(sizeof(struct shape));
    group(s, children, num_children);
    return s;
}

void
group_free(Shape group)
{
    int i;
    for (i = 0; i < group->fields.group.num_children; ++i) {
        shape_free(group->fields.group.children + i);
    }
    free(group->fields.group.children);
    free(group->fields.group.children_xs);
}
