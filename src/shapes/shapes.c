#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libs/linalg/linalg.h"
#include "../pattern/pattern.h"
#include "../material/material.h"
#include "../renderer/ray.h"

#include "group.h"
#include "shapes.h"

Shape
array_of_shapes(size_t num)
{
    return array_of_shapes_realloc(NULL, num);
}

Shape
array_of_shapes_realloc(Shape ptr, size_t num)
{
    return (Shape)realloc(ptr, num * sizeof(struct shape));
}

void
shape_free(Shape s)
{
    if (s != NULL) {
        intersections_free(s->xs);
        s->xs = NULL;
        material_free(s->material);
        s->material = NULL;
        s->bbox_valid = false;
        if (s->type == SHAPE_GROUP) {
            group_free(s);
        }
        // TODO CSG free
    }
}

Intersections
shape_intersect(Shape sh, Ray r, bool stop_after_first_hit)
{
    struct ray transformed_ray;
    Ray r2 = &transformed_ray;

    if (sh->transform_identity) {
        r2 = r;
    } else {
        ray_transform(r, sh->transform_inverse, &transformed_ray);
    }

    Intersections xs = sh->local_intersect(sh, r2, stop_after_first_hit);

    return xs;
}

/*
 * https://docs.unity3d.com/Manual/StandardShaderMaterialParameterNormalMap.html
 * explains how to apply a normal map stored as an RGB image.
 */
void
shape_normal_at(Shape sh, Point world_point, Intersection hit, Vector res)
{
    Point local_point;
    Vector local_normal;
    Vector perturbed;
    Vector world_normal;
    Vector tmp;

    sh->world_to_object(sh, world_point, local_point);
    sh->local_normal_at(sh, local_point, hit, local_normal);
    sh->normal_to_world(sh, local_normal, world_normal);

    if (sh->material->map_bump != NULL) {
        sh->material->map_bump->pattern_at_shape(sh->material->map_bump, sh, world_point, tmp);

        vector_scale(tmp, 2.0);
        perturbed[0] = tmp[0] - 1.0;
        perturbed[1] = tmp[1] - 1.0;
        perturbed[2] = tmp[2] - 1.0;
        perturbed[3] = 0.0;

        //vector_scale(world_normal, 0.5);
        color_accumulate(world_normal, perturbed);
    }

    vector_normalize(world_normal, res);
}

void
shape_normal_to_world(Shape sh, Vector local_normal, Vector res)
{
    Matrix tr;
    Vector un_normal;
    Vector normal;

    matrix_transpose(sh->transform_inverse, tr);

    if (sh->transform_identity) {
        vector_copy(normal, local_normal);
    } else {
        matrix_vector_multiply(tr, local_normal, un_normal);
        vector_normalize(un_normal, normal);
    }

    Vector tmp;
    if (sh->parent != NULL) {
        sh->parent->normal_to_world(sh->parent, normal, tmp);
    } else {
        vector_copy(tmp, normal);
    }
    vector_copy(res, tmp);
}

void
shape_world_to_object(Shape sh, Point pt, Point res)
{
    Point  tmp;
    if (sh->parent != NULL) {
        sh->parent->world_to_object(sh->parent, pt, tmp);
    } else {
        point_copy(tmp, pt);
    }

    if (sh->transform_identity) {
        point_copy(res, pt);
    } else {
        matrix_point_multiply(sh->transform_inverse, tmp, res);
    }
}

void
shape_divide(Shape sh, size_t threshold)
{
    return;
}

bool
shape_includes(Shape a, Shape b)
{
    return a == b;
}

void
shape_set_transform(Shape obj, const Matrix m)
{
    if (obj != NULL) {
        matrix_copy(m, obj->transform);
        matrix_inverse(m, obj->transform_inverse);
        obj->transform_identity = matrix_equal(m, MATRIX_IDENTITY);
    }
}

void
shape_set_material(Shape obj, Material m)
{
    if (obj != NULL) {
        if (obj->material == m) {
            return;
        }
        if (obj->material != NULL) {
            material_free(obj->material);
        }
        obj->material = m;
        if (m != NULL) {
            m->ref_count++;
        }
    }
}

void
shape_set_material_recursive(Shape obj, Material m)
{
    if (obj != NULL) {
        shape_set_material(obj, m);
        if (obj->type == SHAPE_GROUP) {
            int i;
            for (i = 0; i< obj->fields.group.num_children; ++i) {
                shape_set_material_recursive(obj->fields.group.children + i, m);
            }
        }
    }
}

// bounding boxes
void
shape_bounds(Shape sh, Bounding_box *res)
{
    if (!sh->bbox_valid) {
        sh->bbox_valid = true;
        double arr[4] = {-1.0, -1.0, -1.0, 1.0};
        bounding_box_add_array(&(sh->bbox), arr);

        arr[0] = 1.0;
        arr[1] = 1.0;
        arr[2] = 1.0;
        bounding_box_add_array(&(sh->bbox), arr);

        bounding_box_transform(&(sh->bbox), sh->transform, &(sh->bbox_inverse));
    }

    *res = sh->bbox;
}

void
shape_parent_space_bounds(Shape sh, Bounding_box *res)
{
    if (!sh->bbox_valid) {
        sh->bounds(sh, res); // instead of declaring a stack variable
    }

    *res = sh->bbox_inverse;
}

int
shape_to_string(char *buf, size_t n, Shape sh)
{
    const char *type_name;
    switch (sh->type) {
    case SHAPE_CONE:
        type_name = "cone";
        break;
    case SHAPE_CUBE:
        type_name = "cube";
        break;
    case SHAPE_CYLINDER:
        type_name = "cylinder";
        break;
    case SHAPE_PLANE:
        type_name = "plane";
        break;
    case SHAPE_SMOOTH_TRIANGLE:
        type_name = "smooth triangle";
        break;
    case SHAPE_SPHERE:
        type_name = "sphere";
        break;
    case SHAPE_TOROID:
        type_name = "toroid";
        break;
    case SHAPE_TRIANGLE:
        type_name = "triangle";
        break;
    case SHAPE_CSG:
        type_name = "csg";
        break;
    case SHAPE_GROUP:
        type_name = "group";
        break;
    }

    return snprintf(buf, n, "Shape:\n\ttransform: %p\n\tmaterial: %p\n\tparent: %p\n\ttype: %s\n",
        (void *)sh->transform,
        (void *)sh->material,
        (void *)sh->parent,
        type_name);
}


void
shape_recursive_parent_update(Shape sh, Shape parent)
{
    Shape child;
    int i;

    sh->parent = parent;

    if (sh->type == SHAPE_GROUP) {
        for (i = 0, child = sh->fields.group.children;
                i < sh->fields.group.num_children;
                i++, child++) {
            shape_recursive_parent_update(child, sh);
        }
    } else if (sh->type == SHAPE_CSG) {
        shape_recursive_parent_update(sh->fields.csg.left, sh);
        shape_recursive_parent_update(sh->fields.csg.right, sh);
    }
}

void
invalidate_bounding_box(Shape sh)
{
    sh->bbox_valid = false;
    bounding_box(&(sh->bbox));
}

void
shape_recursive_invalidate_bounding_box(Shape sh)
{
    Shape child;
    int i;

    invalidate_bounding_box(sh);

    if (sh->type == SHAPE_GROUP) {
        for (i = 0, child = sh->fields.group.children;
                i < sh->fields.group.num_children;
                i++, child++) {
            shape_recursive_invalidate_bounding_box(child);
        }
    } else if (sh->type == SHAPE_CSG) {
        shape_recursive_invalidate_bounding_box(sh->fields.csg.left);
        shape_recursive_invalidate_bounding_box(sh->fields.csg.right);
    }
}
