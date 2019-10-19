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
        bounding_box_free(s->bbox);
        s->bbox = NULL;
        bounding_box_free(s->bbox_inverse);
        s->bbox_inverse = NULL;
        if (s->type == SHAPE_GROUP) {
            group_free(s);
        }
        // TODO CSG free
    }
}

Intersections
shape_intersect(Shape sh, Ray r)
{
    struct ray transformed_ray;

    ray_transform(r, sh->transform_inverse, &transformed_ray);
    Intersections xs = sh->local_intersect(sh, &transformed_ray);

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
    } else {
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
    
    matrix_vector_multiply(tr, local_normal, un_normal);
    vector_normalize(un_normal, normal);

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

    matrix_point_multiply(sh->transform_inverse, tmp, res);
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
Bounding_box
shape_bounds(Shape sh)
{
    if (sh->bbox == NULL) {
        double arr[4] = {-1.0, -1.0, -1.0, 1.0};

        Bounding_box box = bounding_box_alloc();

        bounding_box_add_array(box, arr);
        arr[0] = 1.0;
        arr[1] = 1.0;
        arr[2] = 1.0;
        bounding_box_add_array(box, arr);

        sh->bbox = box;
        sh->bbox_inverse = bounding_box_transform(box, sh->transform);
    }

    return sh->bbox;
}

Bounding_box
shape_parent_space_bounds(Shape sh)
{
    if (sh->bbox_inverse == NULL) {
        sh->bounds(sh);
    }

    return sh->bbox_inverse;
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
