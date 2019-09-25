#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "shapes.h"
#include "perlin.h"

void
intersection(double t, Shape sh, Intersection x)
{
    x->t = t;
    x->object = sh;
    x->u = -1;
    x->v = -1;
}

void
intersection_with_uv(double t, double u, double v, Shape sh, Intersection x)
{
    x->t = t;
    x->object = sh;
    x->u = u;
    x->v = v;
}

Intersection
intersection_alloc(double t, Shape sh)
{
    Intersection x = (Intersection) malloc(sizeof(struct intersection));
    intersection(t, sh, x);
    return x;
}

Intersection
intersection_uv_alloc(double t, double u, double v, Shape sh)
{
    Intersection x = (Intersection) malloc(sizeof(struct intersection));
    intersection_with_uv(t, u, v, sh, x);
    return x;
}

Intersection
hit(Intersections xs, bool filter_shadow_casters)
{
    int i;
    Intersection x;

    for (i = 0, x = xs->xs; i < xs->num; i++, x++) {
        if (x->t > 0 && (!filter_shadow_casters || (filter_shadow_casters && x->object->material->casts_shadow))) {
            return x;
        }
    }

    return NULL;
}

Intersections
intersections_empty(size_t num)
{
    Intersections xs = (Intersections) malloc(sizeof(struct intersections));
    if (num > 0) {
        xs->xs = (Intersection) malloc(num * sizeof(struct intersection));
    } else {
        xs->xs = NULL;
    }

    xs->array_len = num;
    xs->num = 0;

    return xs;
}

void
intersection_free(Intersection x)
{
    if (x != NULL) {
        free(x);
    }
}

void
intersections_free(Intersections xs)
{
    if (xs != NULL) {
        if (xs->xs != NULL) {
            free(xs->xs);
        }
        free(xs);
    }
}

void
ray_array(double origin[4], double direction[4], Ray ray)
{
    memcpy(ray->origin, origin, 4 * sizeof(double));
    memcpy(ray->direction, direction, 4 * sizeof(double));
}

Ray
ray_array_alloc(double origin[4], double direction[4])
{
    Ray ray = (Ray) malloc(sizeof(struct ray));
    ray_array(origin, direction, ray);
    return ray;
}
Ray
ray_alloc(Point origin, Vector direction)
{
    return ray_array_alloc(origin->arr, direction->arr);
}

void
ray_free(Ray r) {
    if (r != NULL) {
        free(r);
    }
}

void
material(Material m)
{
    m->color[0] = 1.0;
    m->color[1] = 1.0;
    m->color[2] = 1.0;
    m->ambient = 0.1;
    m->diffuse = 0.9;
    m->specular = 0.9;
    m->shininess = 200.0;
    m->reflective = 0.0;
    m->transparency = 0.0;
    m->refractive_index = 1.0;
    m->casts_shadow = true;
    m->pattern = NULL;
}

Material
material_alloc()
{
    Material m = (Material) malloc(sizeof(struct material));
    material(m);
    return m;
}

void
material_set_pattern(Material m, Pattern p)
{
    if (m != NULL) {
        if (m->pattern != NULL) {
            // TODO probably don't need to worry about freeing the pattern ATM
        }
        m->pattern = p;
    }
}

Shape
array_of_shapes(size_t num)
{
    return (Shape)malloc(num * sizeof(struct shape));
}

void
shape_free(Shape s)
{
    if (s != NULL) {
        free(s);
    }
}


void
ray_transform(Ray original, Matrix m, Ray res)
{
    matrix_array_multiply(m, original->origin, res->origin);
    matrix_array_multiply(m, original->direction, res->direction);
}

Ray
ray_transform_alloc(Ray original, Matrix m)
{
    Ray res = ray_array_alloc(original->origin, original->direction);
    ray_transform(original, m, res);
    return res;
}

int
ray_to_string(char *buf, size_t n, Ray r)
{
    return snprintf(buf, n, "Point: [%f %f %f] Vector: [%f %f %f]",
                    r->origin[0], r->origin[1], r->origin[2],
                    r->direction[0], r->direction[1], r->direction[2]);
}

Intersections
shape_intersect(Shape sh, Ray r)
{
    Matrix m = sh->transform_inverse;
    struct ray transformed_ray;

    ray_transform(r, m, &transformed_ray);
    Intersections xs = sh->local_intersect(sh, &transformed_ray);

    return xs;
}

Vector
shape_normal_at(Shape sh, Point world_point, Intersection hit)
{
    Point local_point = sh->world_to_object(sh, world_point);
    Vector local_normal = sh->local_normal_at(sh, local_point, hit);
    Vector world_normal = sh->normal_to_world(sh, local_normal);
    Vector n = vector_normalize_alloc(world_normal);

    vector_free(world_normal);
    vector_free(local_normal);
    point_free(local_point);

    return n;
}

Vector
shape_normal_to_world(Shape sh, Vector local_normal)
{
    Matrix inv = sh->transform_inverse;
    Matrix tr = matrix_transpose_alloc(inv);
    
    Vector un_normal = matrix_vector_multiply_alloc(tr, local_normal);
    
    Vector normal = vector_normalize_alloc(un_normal);

    if (sh->parent != NULL) {
        Vector new_normal = sh->parent->normal_to_world(sh->parent, normal);
        vector_free(normal);
        normal = new_normal;
    }

    vector_free(un_normal);
    matrix_free(tr);

    return normal;
}

Point
shape_world_to_object(Shape sh, Point pt)
{
    Point pt_new = pt;
    if (sh->parent != NULL) {
        pt_new = sh->parent->world_to_object(sh->parent, pt);
    }

    Matrix m = sh->transform_inverse;
    Point p = matrix_point_multiply_alloc(m, pt_new);

    if (sh->parent != NULL) {
        point_free(pt_new);
    }

    return p;
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
shape_set_transform(Shape obj, Matrix m)
{
    if (obj != NULL) {
        if (obj->transform_inverse != NULL) {
            matrix_free(obj->transform_inverse);
        }
        if (obj->transform != NULL) {
            matrix_free(obj->transform);
        }
        obj->transform = m;
        obj->transform_inverse = matrix_inverse_alloc(m);
    }
}

void
shape_set_material(Shape obj, Material m)
{
    if (obj != NULL) {
        obj->material = m;
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

int
intersection_to_string(char *buf, size_t n, Intersection x)
{
    int used = snprintf(buf, n, "Intersection:\n\tt: %f\n\tu: %f\n\tv: %f\n",
                        x->t,
                        x->u,
                        x->v);
    return used + shape_to_string(buf + used, n - used, x->object);
}

/*
 *
 * Patterns
 *
 */
Color
base_pattern_at_shape(Pattern p, Shape s, Point pt)
{
    Matrix inv = p->transform_inverse;

    Point object_point = s->world_to_object(s, pt);
    Point pattern_point = matrix_point_multiply_alloc(inv, object_point);
    Color c = p->pattern_at(p, pattern_point);

    point_free(pattern_point);
    point_free(object_point);

    return c;
}

Color
blended_pattern_at_shape(Pattern p, Shape s, Point pt)
{
    Color c1 = p->fields.blended.pattern1->pattern_at_shape(p->fields.blended.pattern1, s, pt);
    Color c2 = p->fields.blended.pattern2->pattern_at_shape(p->fields.blended.pattern2, s, pt);

    Color c3 = color((c1->arr[0] + c2->arr[0]) / 2,
                     (c1->arr[1] + c2->arr[1]) / 2,
                     (c1->arr[2] + c2->arr[2]) / 2);

    color_free(c1);
    color_free(c2);

    return c3;
}

Color
nested_pattern_at_shape(Pattern p, Shape s, Point pt)
{
    Color c1 = p->fields.nested.pattern2->pattern_at_shape(p->fields.nested.pattern2, s, pt);
    Color c2 = p->fields.nested.pattern3->pattern_at_shape(p->fields.nested.pattern3, s, pt);

    switch(p->fields.nested.pattern1->type) {
    case CHECKER_PATTERN:
    case GRADIENT_PATTERN:
    case RADIAL_GRADIENT_PATTERN:
    case RING_PATTERN:
    case STRIPE_PATTERN:
        memcpy(p->fields.nested.pattern1->fields.concrete.a, c1->arr, 3 * sizeof(double));
        memcpy(p->fields.nested.pattern1->fields.concrete.b, c2->arr, 3 * sizeof(double));
        break;
    // TODO I can't think of anything to do for these other patterns...
    case UV_ALIGN_CHECKER_PATTERN:
    case UV_CHECKER_PATTERN:
    case UV_TEXTURE_PATTERN:
    case BLENDED_PATTERN:
    case NESTED_PATTERN:
    case PERTURBED_PATTERN:
    case CUBE_MAP_PATTERN:
    case CYLINDER_MAP_PATTERN:
    case TEXTURE_MAP_PATTERN:
    default:
        break;
    }

    color_free(c1);
    color_free(c2);

    return p->fields.nested.pattern1->pattern_at_shape(p->fields.nested.pattern1, s, pt);
}

Color
perturbed_pattern_at_shape(Pattern p, Shape s, Point pt)
{
    double x = pt->arr[0] / 10.0;
    double y = pt->arr[1] / 10.0;
    double z = pt->arr[2] / 10.0;

    double new_x = pt->arr[0] + p->fields.perturbed.scale_factor *
                   pnoise3d(x, y, z,
                            p->fields.perturbed.persistence,
                            p->fields.perturbed.frequency,
                            p->fields.perturbed.octaves,
                            p->fields.perturbed.seed);

    z += 1.0;
    double new_y = pt->arr[1] + p->fields.perturbed.scale_factor *
                   pnoise3d(x, y, z,
                            p->fields.perturbed.persistence,
                            p->fields.perturbed.frequency,
                            p->fields.perturbed.octaves,
                            p->fields.perturbed.seed);

    z += 1.0;
    double new_z = pt->arr[2] + p->fields.perturbed.scale_factor *
                   pnoise3d(x, y, z,
                            p->fields.perturbed.persistence,
                            p->fields.perturbed.frequency,
                            p->fields.perturbed.octaves,
                            p->fields.perturbed.seed);

    Point perturbed = point(new_x, new_y, new_z);
    Color c = p->fields.perturbed.pattern1->pattern_at_shape(p->fields.perturbed.pattern1, s, perturbed);

    point_free(perturbed);

    return c;
}

Color
base_pattern_at(Pattern p, Point pt)
{
    printf("calling base_pattern_at which should only happen for tests.\n");
    return color(pt->arr[0], pt->arr[1], pt->arr[2]);
}

Color
checker_pattern_at(Pattern p, Point pt)
{
    int s = (int)floor(pt->arr[0]) + (int)floor(pt->arr[1]) + (int)floor(pt->arr[2]);
    double *arr;

    if (s % 2 == 0) {
        arr = p->fields.concrete.a;
    } else {
        arr = p->fields.concrete.b;
    }
    return color(arr[0], arr[1], arr[2]);
}

Color
gradient_pattern_at(Pattern p, Point pt)
{
    double distance[3] = { p->fields.concrete.b[0] - p->fields.concrete.a[0],
                           p->fields.concrete.b[1] - p->fields.concrete.a[1],
                           p->fields.concrete.b[2] - p->fields.concrete.a[2]};
    double fraction = pt->arr[0] - floor(pt->arr[0]);

    return color(p->fields.concrete.a[0] + distance[0] * fraction,
                 p->fields.concrete.a[1] + distance[1] * fraction,
                 p->fields.concrete.a[2] + distance[2] * fraction);
}

Color
radial_gradient_pattern_at(Pattern p, Point pt)
{
    double distance[3] = { p->fields.concrete.b[0] - p->fields.concrete.a[0],
                           p->fields.concrete.b[1] - p->fields.concrete.a[1],
                           p->fields.concrete.b[2] - p->fields.concrete.a[2]};

    double magnitude = sqrt(pt->arr[0] * pt->arr[0] + pt->arr[2] * pt->arr[2]);

    double fraction = magnitude - floor(magnitude);

    return color(p->fields.concrete.a[0] + distance[0] * fraction,
                 p->fields.concrete.a[1] + distance[1] * fraction,
                 p->fields.concrete.a[2] + distance[2] * fraction);
}

Color
ring_pattern_at(Pattern p, Point pt)
{
    int s = (int)floor(sqrt(pt->arr[0] * pt->arr[0] + pt->arr[2] * pt->arr[2]));
    double *arr;

    if (s % 2 == 0) {
        arr = p->fields.concrete.a;
    } else {
        arr = p->fields.concrete.b;
    }
    return color(arr[0], arr[1], arr[2]);
}

Color
stripe_pattern_at(Pattern p, Point pt)
{
    int s = (int)floor(pt->arr[0]);
    double *arr;

    if (s % 2 == 0) {
        arr = p->fields.concrete.a;
    } else {
        arr = p->fields.concrete.b;
    }
    return color(arr[0], arr[1], arr[2]);
}


UVMapReturnType
uv_map_return_type_alloc(size_t face, double u, double v)
{
    UVMapReturnType rv = (UVMapReturnType) malloc(sizeof(struct face_uv_retval));
    rv->face = face;
    rv->u = u;
    rv->v = v;

    return rv;
}

void
uv_map_return_type_free(UVMapReturnType p)
{
    if (p != NULL) {
        free(p);
    }
}


Color
texture_map_pattern_at(Pattern p, Point pt)
{
    UVMapReturnType face_u_v = p->uv_map(pt);
    Pattern face = p->fields.uv_map.uv_faces + face_u_v->face;

    Color c = face->uv_pattern_at(face, face_u_v->u, face_u_v->v);

    uv_map_return_type_free(face_u_v);

    return c;
}

Color
base_uv_pattern_at(Pattern p, double u, double v)
{
    printf("calling base_uv_pattern_at which should never happen.\n");
    return color(u, v, 0);
}

Color
uv_align_check_uv_pattern_at(Pattern p, double u, double v)
{
    // remember v=0 at the bottom, v=1 at the top
    double *arr = p->fields.uv_align_check.main;

    if (v > 0.8) {
        if (u < 0.2) {
            arr = p->fields.uv_align_check.ul;
        } else if (u > 0.8) {
            arr = p->fields.uv_align_check.ur;
        }
    } else if (v < 0.2) {
        if (u < 0.2) {
            arr = p->fields.uv_align_check.bl;
        } else if (u > 0.8) {
            arr = p->fields.uv_align_check.br;
        }
    }

    return color(arr[0], arr[1], arr[2]);
}

Color
uv_check_uv_pattern_at(Pattern p, double u, double v)
{
    double *arr;
    int u2 = (int)floor(u * (double)p->fields.uv_check.width),
        v2 = (int)floor(v * (double)p->fields.uv_check.height);

    if ((u2 + v2) % 2 == 0) {
        arr = p->fields.uv_check.a;
    } else {
        arr = p->fields.uv_check.b;
    }

    return color(arr[0], arr[1], arr[2]);
}

Color
uv_texture_uv_pattern_at(Pattern p, double u, double v)
{
    // remember v=0 at the bottom, v=1 at the top
    v = 1 - v;
    int y = (int)round(u * (double)(p->fields.uv_texture.canvas->width - 1));
    int x = (int)round(v * (double)(p->fields.uv_texture.canvas->height - 1));

    Color c = color(0,0,0);
    canvas_pixel_at(p->fields.uv_texture.canvas, y, x, c);

    return c;
}

UVMapReturnType
base_uv_map(Point pt)
{
    printf("calling base_uv_pattern_at which should never happen.\n");
    return uv_map_return_type_alloc(0, pt->arr[0], pt->arr[1]);
}

UVMapReturnType
cube_uv_map(Point pt)
{
    double abs_x = fabs(pt->arr[0]);
    double abs_y = fabs(pt->arr[1]);
    double abs_z = fabs(pt->arr[2]);
    double coord = fmax(fmax(abs_x, abs_y), abs_z);

    size_t face = equal(coord, pt->arr[0])
        ? 0
        : equal(coord, -pt->arr[0])
            ? 1
            : equal(coord, pt->arr[1])
                ? 2
                : equal(coord, -pt->arr[1])
                    ? 3
                    : equal(coord, pt->arr[2])
                        ? 4
                        : 5;

    UVMapReturnType retval = uv_map_return_type_alloc(face, 0, 0);

    switch (face) {
    case 0: // right
        retval->u = fmod((1.0 - pt->arr[2]), 2.0) / 2.0;
        retval->v = fmod((pt->arr[1] + 1.0), 2.0) / 2.0;
        break;
    case 1: // left
        retval->u = fmod((pt->arr[2] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt->arr[1] + 1.0), 2.0) / 2.0;
        break;
    case 2: // up
        retval->u = fmod((pt->arr[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((1.0 - pt->arr[2]), 2.0) / 2.0;
        break;
    case 3: // down
        retval->u = fmod((pt->arr[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt->arr[2] + 1.0), 2.0) / 2.0;
        break;
    case 4: // front
        retval->u = fmod((pt->arr[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt->arr[1] + 1.0), 2.0) / 2.0;
        break;
    case 5: // back
        retval->u = fmod((1.0 - pt->arr[0]), 2.0) / 2.0;
        retval->v = fmod((pt->arr[1] + 1.0), 2.0) / 2.0;
        break;
    }

    return retval;
}

UVMapReturnType
cylinder_uv_map(Point pt)
{
    //double radius = sqrt(pt->arr[0] * pt->arr[0] + pt->arr[2] + pt->arr[2]);
    //double abs_y = fabs(pt->arr[1]);
    //double coord = fmax(radius, abs_y);
    double theta;
    double raw_u;

    // TODO for now don't bother with special mapping for the caps
    size_t face = 0; /*equal(coord, radius)
        ? 0
        : equal(coord, pt->arr[1])
            ? 1
            : 2; */

    UVMapReturnType retval = uv_map_return_type_alloc(face, 0, 0);

    switch (face) {
    case 0: // cylinder body
        theta = atan2(pt->arr[0], pt->arr[2]);
        raw_u = theta / (2.0 * M_PI);
        retval->u = 1.0 - (raw_u + 0.5);
        retval->v = fmod(pt->arr[1], 1.0);
        break;
    case 1: // top cap
        retval->u = fmod((pt->arr[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((1.0 - pt->arr[2]), 2.0) / 2.0;
        break;
    case 2: // bottom cap
        retval->u = fmod((pt->arr[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt->arr[2] + 1.0), 2.0) / 2.0;
        break;
    }

    return retval;
}

UVMapReturnType
plane_uv_map(Point pt)
{
    double u = fmod(pt->arr[0], 1.0);
    double v = fmod(pt->arr[2], 1.0);
    if (u < 0) {
        u += 1.0;
    }
    if (v < 0) {
        v += 1.0;
    }

    return uv_map_return_type_alloc(0, u, v);
}

UVMapReturnType
sphere_uv_map(Point pt)
{
    double theta = atan2(pt->arr[0], pt->arr[2]);
    Vector vec = vector(pt->arr[0], pt->arr[1], pt->arr[2]);
    double radius = vector_magnitude(vec);
    double phi = acos(pt->arr[1] / radius);
    double raw_u = theta / (2 * M_PI);
    double u = 1 - (raw_u + 0.5);
    double v = 1 - phi / M_PI;

    vector_free(vec);

    return uv_map_return_type_alloc(0, u, v);
}


void
pattern_set_transform(Pattern obj, Matrix m)
{
    if (obj != NULL) {
        if (obj->transform_inverse != NULL) {
            matrix_free(obj->transform_inverse);
        }
        if (obj->transform != NULL) {
            matrix_free(obj->transform);
        }
        obj->transform = m;
        obj->transform_inverse = matrix_inverse_alloc(m);
    }
}

void
default_pattern_constructor(Pattern res)
{
    res->transform = NULL;
    res->transform_inverse = NULL;

    pattern_set_transform(res, matrix_identity_alloc());

    res->pattern_at_shape = base_pattern_at_shape;
    res->pattern_at = base_pattern_at;
    res->uv_pattern_at = base_uv_pattern_at;
    res->uv_map = base_uv_map;
}

void
checker_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = CHECKER_PATTERN;
    memcpy(res->fields.concrete.a, a->arr, 3 * sizeof(double));
    memcpy(res->fields.concrete.b, b->arr, 3 * sizeof(double));

    res->pattern_at = checker_pattern_at;
}

void
gradient_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = GRADIENT_PATTERN;
    memcpy(res->fields.concrete.a, a->arr, 3 * sizeof(double));
    memcpy(res->fields.concrete.b, b->arr, 3 * sizeof(double));

    res->pattern_at = gradient_pattern_at;
}

void
radial_gradient_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = RADIAL_GRADIENT_PATTERN;
    memcpy(res->fields.concrete.a, a->arr, 3 * sizeof(double));
    memcpy(res->fields.concrete.b, b->arr, 3 * sizeof(double));

    res->pattern_at = radial_gradient_pattern_at;
}

void
ring_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = RING_PATTERN;
    memcpy(res->fields.concrete.a, a->arr, 3 * sizeof(double));
    memcpy(res->fields.concrete.b, b->arr, 3 * sizeof(double));

    res->pattern_at = ring_pattern_at;
}

void
stripe_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = STRIPE_PATTERN;
    memcpy(res->fields.concrete.a, a->arr, 3 * sizeof(double));
    memcpy(res->fields.concrete.b, b->arr, 3 * sizeof(double));

    res->pattern_at = stripe_pattern_at;
}

void
uv_align_check_pattern(Color main, Color ul, Color ur, Color bl, Color br, Pattern res)
{
    default_pattern_constructor(res);
    res->type = UV_ALIGN_CHECKER_PATTERN;

    memcpy(res->fields.uv_align_check.main, main->arr, 3 * sizeof(double));
    memcpy(res->fields.uv_align_check.ul, ul->arr, 3 * sizeof(double));
    memcpy(res->fields.uv_align_check.ur, ur->arr, 3 * sizeof(double));
    memcpy(res->fields.uv_align_check.bl, bl->arr, 3 * sizeof(double));
    memcpy(res->fields.uv_align_check.br, br->arr, 3 * sizeof(double));

    res->uv_pattern_at = uv_align_check_uv_pattern_at;
}

void
uv_check_pattern(Color a, Color b, size_t width, size_t height, Pattern res)
{
    default_pattern_constructor(res);
    res->type = UV_CHECKER_PATTERN;

    memcpy(res->fields.uv_check.a, a->arr, 3 * sizeof(double));
    memcpy(res->fields.uv_check.b, b->arr, 3 * sizeof(double));
    res->fields.uv_check.width = width;
    res->fields.uv_check.height = height;

    res->uv_pattern_at = uv_check_uv_pattern_at;
}

void
uv_texture_pattern(Canvas canvas, Pattern res)
{
    default_pattern_constructor(res);
    res->type = UV_TEXTURE_PATTERN;

    res->fields.uv_texture.canvas = canvas;

    res->uv_pattern_at = uv_texture_uv_pattern_at;
}

void
blended_pattern(Pattern p1, Pattern p2, Pattern res)
{
    default_pattern_constructor(res);
    res->type = BLENDED_PATTERN;

    res->fields.blended.pattern1 = p1;
    res->fields.blended.pattern2 = p2;

    res->pattern_at_shape = blended_pattern_at_shape;
}

void
nested_pattern(Pattern p1, Pattern p2, Pattern p3, Pattern res)
{
    default_pattern_constructor(res);
    res->type = NESTED_PATTERN;

    res->fields.nested.pattern1 = p1;
    res->fields.nested.pattern2 = p2;
    res->fields.nested.pattern3 = p3;

    res->pattern_at_shape = nested_pattern_at_shape;
}

void
perturbed_pattern(Pattern p1, double frequency, double scale_factor, double persistence, size_t octaves, int seed, Pattern res)
{
    default_pattern_constructor(res);
    res->type = PERTURBED_PATTERN;

    res->fields.perturbed.pattern1 = p1;
    res->fields.perturbed.frequency = frequency;
    res->fields.perturbed.scale_factor = scale_factor;
    res->fields.perturbed.persistence = persistence;
    res->fields.perturbed.octaves = octaves;
    res->fields.perturbed.seed = seed;

    res->pattern_at_shape = perturbed_pattern_at_shape;
}

Pattern
array_of_patterns(size_t num)
{
    return (Pattern) malloc(num * sizeof(struct pattern));
}

void
texture_map_pattern(Pattern faces, enum uv_map_type type, Pattern res)
{
    default_pattern_constructor(res);
    res->type = TEXTURE_MAP_PATTERN;

    res->fields.uv_map.type = type;
    res->fields.uv_map.uv_faces = faces;
    res->pattern_at = texture_map_pattern_at;

    switch (type) {
    case CUBE_UV_MAP:
        res->uv_map = cube_uv_map;
        break;
    case CYLINDER_UV_MAP:
        res->uv_map = cylinder_uv_map;
        break;
    case PLANE_UV_MAP:
        res->uv_map = plane_uv_map;
        break;
    case SPHERE_UV_MAP:
        res->uv_map = sphere_uv_map;
        break;
    default:
        // already set res->uv_map to an error function in the default constructor
        break;
    }
}

Pattern
checker_pattern_alloc(Color a, Color b)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    checker_pattern(a, b, p);
    return p;
}

Pattern
gradient_pattern_alloc(Color a, Color b)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    gradient_pattern(a, b, p);
    return p;
}

Pattern
radial_gradient_pattern_alloc(Color a, Color b)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    radial_gradient_pattern(a, b, p);
    return p;
}

Pattern ring_pattern_alloc(Color a, Color b)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    ring_pattern(a, b, p);
    return p;
}

Pattern
stripe_pattern_alloc(Color a, Color b)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    stripe_pattern(a, b, p);
    return p;
}

Pattern
uv_align_check_pattern_alloc(Color main, Color ul, Color ur, Color bl, Color br)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    uv_align_check_pattern(main, ul, ur, bl, br, p);
    return p;
}

Pattern
uv_check_pattern_alloc(Color a, Color b, size_t width, size_t height)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    uv_check_pattern(a, b, width, height, p);
    return p;
}

Pattern
uv_texture_pattern_alloc(Canvas canvas)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    uv_texture_pattern(canvas, p);
    return p;
}

Pattern
blended_pattern_alloc(Pattern p1, Pattern p2)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    blended_pattern(p1, p2, p);
    return p;
}

Pattern
nested_pattern_alloc(Pattern p1, Pattern p2, Pattern p3)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    nested_pattern(p1, p2, p3, p);
    return p;
}

Pattern
perturbed_pattern_alloc(Pattern p1, double frequency, double scale_factor, double persistence, size_t octaves, int seed)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    perturbed_pattern(p1, frequency, scale_factor, persistence, octaves, seed, p);
    return p;
}

Pattern
texture_map_pattern_alloc(Pattern faces /* number should be determined by uv_map */,  enum uv_map_type type)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    texture_map_pattern(faces, type, p);
    return p;
}


