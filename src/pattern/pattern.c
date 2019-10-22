#include <math.h>
#include <string.h>

#include "../libs/perlin/perlin.h"
#include "../color/color.h"

#include "pattern.h"

void
base_pattern_at_shape(Pattern p, Shape s, Point pt, Color res)
{
    Point object_point;
    Point pattern_point;
    Color c;

    s->world_to_object(s, pt, object_point);

    if (p->transform_identity) {
        point_copy(pattern_point, object_point);
    } else {
        matrix_point_multiply(p->transform_inverse, object_point, pattern_point);
    }

    p->pattern_at(p, s, pattern_point, c);

    color_copy(res, c);
    res[3] = 0.0;
}

void
blended_pattern_at_shape(Pattern p, Shape s, Point pt, Color res)
{
    Color c1;
    Color c2;
    p->fields.blended.pattern1->pattern_at_shape(p->fields.blended.pattern1, s, pt, c1);
    p->fields.blended.pattern2->pattern_at_shape(p->fields.blended.pattern2, s, pt, c2);

    color_average(c1, c2, res);
}

// TODO not threadsafe
void
nested_pattern_at_shape(Pattern p, Shape s, Point pt, Color res)
{
    Color c1;
    Color c2;
    Color c3;
    p->fields.nested.pattern2->pattern_at_shape(p->fields.nested.pattern2, s, pt, c1);
    p->fields.nested.pattern3->pattern_at_shape(p->fields.nested.pattern3, s, pt, c2);

    switch(p->fields.nested.pattern1->type) {
    case CHECKER_PATTERN:
    case GRADIENT_PATTERN:
    case RADIAL_GRADIENT_PATTERN:
    case RING_PATTERN:
    case STRIPE_PATTERN:
        color_copy(p->fields.nested.pattern1->fields.concrete.a, c1);
        color_copy(p->fields.nested.pattern1->fields.concrete.b, c2);
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

    p->fields.nested.pattern1->pattern_at_shape(p->fields.nested.pattern1, s, pt, c3);
    color_copy(res, c3);
}

void
perturbed_pattern_at_shape(Pattern p, Shape s, Point pt, Color res)
{
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];

    double new_x = pt[0] + p->fields.perturbed.scale_factor *
                   pnoise3d(x, y, z,
                            p->fields.perturbed.persistence,
                            p->fields.perturbed.frequency,
                            p->fields.perturbed.octaves,
                            p->fields.perturbed.seed);

    if (z < 0) { z -= 1.0; } else { z += 1.0; }
    double new_y = pt[1] + p->fields.perturbed.scale_factor *
                   pnoise3d(x, y, z,
                            p->fields.perturbed.persistence,
                            p->fields.perturbed.frequency,
                            p->fields.perturbed.octaves,
                            p->fields.perturbed.seed);

    if (z < 0) { z -= 1.0; } else { z += 1.0; }
    double new_z = pt[2] + p->fields.perturbed.scale_factor *
                   pnoise3d(x, y, z,
                            p->fields.perturbed.persistence,
                            p->fields.perturbed.frequency,
                            p->fields.perturbed.octaves,
                            p->fields.perturbed.seed);

    Point  perturbed;
    perturbed[0] = new_x;
    perturbed[1] = new_y;
    perturbed[2] = new_z;

    Color c;
    p->fields.perturbed.pattern1->pattern_at_shape(p->fields.perturbed.pattern1, s, perturbed, c);
    color_copy(res, c);
}

void
base_pattern_at(Pattern p, Shape s, Point pt, Color res)
{
    printf("calling base_pattern_at which should only happen for tests.\n");
    color_copy(res, pt);
}

void
checker_pattern_at(Pattern p, Shape s, Point pt, Color res)
{
    int t = (int)floor(pt[0]) + (int)floor(pt[1]) + (int)floor(pt[2]);
    double *arr;

    if (t % 2 == 0) {
        arr = p->fields.concrete.a;
    } else {
        arr = p->fields.concrete.b;
    }

    color_copy(res, arr);
}

void
gradient_pattern_at(Pattern p, Shape s, Point pt, Color res)
{
    double distance[3] = { p->fields.concrete.b[0] - p->fields.concrete.a[0],
                           p->fields.concrete.b[1] - p->fields.concrete.a[1],
                           p->fields.concrete.b[2] - p->fields.concrete.a[2]};
    double fraction = pt[0] - floor(pt[0]);

    res[0] = p->fields.concrete.a[0] + distance[0] * fraction;
    res[1] = p->fields.concrete.a[1] + distance[1] * fraction;
    res[2] = p->fields.concrete.a[2] + distance[2] * fraction;
}

void
radial_gradient_pattern_at(Pattern p, Shape s, Point pt, Color res)
{
    double distance[3] = { p->fields.concrete.b[0] - p->fields.concrete.a[0],
                           p->fields.concrete.b[1] - p->fields.concrete.a[1],
                           p->fields.concrete.b[2] - p->fields.concrete.a[2]};

    double magnitude = sqrt(pt[0] * pt[0] + pt[2] * pt[2]);

    double fraction = magnitude - floor(magnitude);

    res[0] = p->fields.concrete.a[0] + distance[0] * fraction;
    res[1] = p->fields.concrete.a[1] + distance[1] * fraction;
    res[2] = p->fields.concrete.a[2] + distance[2] * fraction;
}

void
ring_pattern_at(Pattern p, Shape s, Point pt, Color res)
{
    int t = (int)floor(sqrt(pt[0] * pt[0] + pt[2] * pt[2]));
    double *arr;

    if (t % 2 == 0) {
        arr = p->fields.concrete.a;
    } else {
        arr = p->fields.concrete.b;
    }
    color_copy(res, arr);
}

void
stripe_pattern_at(Pattern p, Shape s, Point pt, Color res)
{
    int t = (int)floor(pt[0]);
    double *arr;

    if (t % 2 == 0) {
        arr = p->fields.concrete.a;
    } else {
        arr = p->fields.concrete.b;
    }
    color_copy(res, arr);
}

void
texture_map_pattern_at(Pattern p, Shape s, Point pt, Color res)
{
    UVMapReturnType face_u_v;
    Point uv_pattern_point;
    Color c;

    // call this twice. First to get the face, second with the face's transform applied
    p->uv_map(s, pt, &face_u_v);
    Pattern face = p->fields.uv_map.uv_faces + face_u_v.face;

    if (face->transform_identity) {
        point_copy(uv_pattern_point, pt);
    } else {
        matrix_point_multiply(face->transform_inverse, pt, uv_pattern_point);
    }
    p->uv_map(s, uv_pattern_point, &face_u_v);

    face->uv_pattern_at(face, face_u_v.u, face_u_v.v, c);
    color_copy(res, c);
}

void
base_uv_pattern_at(Pattern p, double u, double v, Color res)
{
    printf("calling base_uv_pattern_at which should never happen.\n");
    res[0] = u;
    res[1] = v;
    res[2] = 0;
}

void
uv_align_check_uv_pattern_at(Pattern p, double u, double v, Color res)
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

    color_copy(res, arr);
}

void
uv_check_uv_pattern_at(Pattern p, double u, double v, Color res)
{
    double *arr;
    int u2 = (int)floor(u * (double)p->fields.uv_check.width),
        v2 = (int)floor(v * (double)p->fields.uv_check.height);

    if ((u2 + v2) % 2 == 0) {
        arr = p->fields.uv_check.a;
    } else {
        arr = p->fields.uv_check.b;
    }

    color_copy(res, arr);
}

void gradient_pattern_at(Pattern p, Shape s, Point pt, Color res);

void
uv_gradient_uv_pattern_at(Pattern p, double u, double v, Color res)
{
    Point point = {u, v, 0.0, 1.0};
    gradient_pattern_at(p, NULL, point, res);
}

void radial_gradient_pattern_at(Pattern p, Shape s, Point pt, Color res);

void
uv_radial_gradient_uv_pattern_at(Pattern p, double u, double v, Color res)
{
    Point point = {u, v, 0.0, 1.0};
    radial_gradient_pattern_at(p, NULL, point, res);
}


void
uv_texture_uv_pattern_at(Pattern p, double u, double v, Color res)
{
    // remember v=0 at the bottom, v=1 at the top
    v = 1 - v;
    size_t y = (size_t)round(u * (double)(p->fields.uv_texture.canvas->width - 1));
    size_t x = (size_t)round(v * (double)(p->fields.uv_texture.canvas->height - 1));

    Color c;
    canvas_pixel_at(p->fields.uv_texture.canvas, y, x, c);

    color_copy(res, c);
}

void
base_uv_map(Shape s, Point pt, UVMapReturnType *retval)
{
    printf("calling base_uv_pattern_at which should never happen.\n");
    retval->face = 0;
    retval->u = pt[0];
    retval->v = pt[1];
}

void
cube_uv_map(Shape s, Point pt, UVMapReturnType *retval)
{
    double abs_x = fabs(pt[0]);
    double abs_y = fabs(pt[1]);
    double abs_z = fabs(pt[2]);
    double coord = fmax(fmax(abs_x, abs_y), abs_z);

    size_t face = equal(coord, pt[0])
        ? 0
        : equal(coord, -pt[0])
            ? 1
            : equal(coord, pt[1])
                ? 2
                : equal(coord, -pt[1])
                    ? 3
                    : equal(coord, pt[2])
                        ? 4
                        : 5;

    retval->face = face;

    switch (face) {
    case 0: // right
        retval->u = fmod((1.0 - pt[2]), 2.0) / 2.0;
        retval->v = fmod((pt[1] + 1.0), 2.0) / 2.0;
        break;
    case 1: // left
        retval->u = fmod((pt[2] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt[1] + 1.0), 2.0) / 2.0;
        break;
    case 2: // up
        retval->u = fmod((pt[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((1.0 - pt[2]), 2.0) / 2.0;
        break;
    case 3: // down
        retval->u = fmod((pt[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt[2] + 1.0), 2.0) / 2.0;
        break;
    case 4: // front
        retval->u = fmod((pt[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt[1] + 1.0), 2.0) / 2.0;
        break;
    case 5: // back
        retval->u = fmod((1.0 - pt[0]), 2.0) / 2.0;
        retval->v = fmod((pt[1] + 1.0), 2.0) / 2.0;
        break;
    }
}

void
cylinder_uv_map(Shape s, Point pt, UVMapReturnType *retval)
{
    double theta;
    double raw_u;

    size_t face = 
        (s->fields.cylinder.maximum - EPSILON) <= pt[1]
        ? 1
        : (s->fields.cylinder.minimum + EPSILON) >= pt[1]
            ? 2
            : 0;

    retval->face = face;

    switch (face) {
    case 0: // cylinder body
        theta = atan2(pt[0], pt[2]);
        raw_u = theta / (2.0 * M_PI);
        retval->u = 1.0 - (raw_u + 0.5);
        retval->v = fmod(pt[1], 1.0);
        break;
    case 1: // top cap
        retval->u = fmod((pt[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((1.0 - pt[2]), 2.0) / 2.0;
        break;
    case 2: // bottom cap
        retval->u = fmod((pt[0] + 1.0), 2.0) / 2.0;
        retval->v = fmod((pt[2] + 1.0), 2.0) / 2.0;
        break;
    }
}

void
triangle_uv_map(Shape s, Point pt, UVMapReturnType *retval)
{
    retval->face = 0;

    Vector v1, v2, v3;

    // https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    // Compute barycentric coordinates (u, v, w) for point pt with respect to triangle (a, b, c)
    vector_from_points(pt, s->fields.triangle.p1, v2);
    
    double d00 = vector_dot(s->fields.triangle.e1, s->fields.triangle.e1);
    double d01 = vector_dot(s->fields.triangle.e1, s->fields.triangle.e2);
    double d11 = vector_dot(s->fields.triangle.e2, s->fields.triangle.e2);
    double d20 = vector_dot(v2, s->fields.triangle.e1);
    double d21 = vector_dot(v2, s->fields.triangle.e2);
    double denom = 1.0 / (d00 * d11 - d01 * d01);
    double u, v, w;
    v = fmod((d11 * d20 - d01 * d21) * denom, 1.0);
    w = fmod((d00 * d21 - d01 * d20) * denom, 1.0);
    u = 1.0 - v - w;

    if (s->fields.triangle.use_textures) {
        vector_copy(v1, s->fields.triangle.t1);
        vector_copy(v2, s->fields.triangle.t2);
        vector_copy(v3, s->fields.triangle.t3);

        vector_scale(v1, u);
        vector_scale(v2, v);
        vector_scale(v3, 1.0 - u - v);

        v1[0] += v2[0] + v3[0];
        v1[1] += v2[1] + v3[1];
        v1[2] += v2[2] + v3[2];

        retval->u = fmod(v1[0], 1.0);
        retval->v = fmod(v1[1], 1.0);
    } else {
        retval->u = u;
        retval->v = v;
    }

    if (retval->u < 0) {
        retval->u += 1.0;
    }
    if (retval->v < 0) {
        retval->v += 1.0;
    }
}

void
plane_uv_map(Shape s, Point pt, UVMapReturnType *retval)
{
    double u = fmod(pt[0], 1.0);
    double v = fmod(pt[2], 1.0);
    if (u < 0) {
        u += 1.0;
    }
    if (v < 0) {
        v += 1.0;
    }

    retval->face = 0;
    retval->u = u;
    retval->v = v;
}

void
sphere_uv_map(Shape s, Point pt, UVMapReturnType *retval)
{
    double theta = atan2(pt[0], pt[2]);
    Vector vec;
    vector_copy(vec, pt);
    vec[3] = 0.0;
    double radius = vector_magnitude(vec);
    double phi = acos(pt[1] / radius);
    double raw_u = theta / (2 * M_PI);
    double u = 1 - (raw_u + 0.5);
    double v = 1 - phi / M_PI;

    retval->face = 0;
    retval->u = u;
    retval->v = v;
}

void
toroid_uv_map(Shape s, Point pt, UVMapReturnType *retval)
{
    double u = (1.0 - (atan2(pt[2], pt[0]) + M_PI) / (2 * M_PI));
    double len = sqrt(pt[0] * pt[0] + pt[2] * pt[2]);
    double x = len - s->fields.toroid.r1;
    double v = (atan2(pt[1], x) + M_PI) / (2 * M_PI);

    retval->face = 0;
    retval->u = u;
    retval->v = v;
}


void
pattern_set_transform(Pattern obj, const Matrix m)
{
    if (obj != NULL) {
        matrix_copy(m, obj->transform);
        matrix_inverse(m, obj->transform_inverse);
        obj->transform_identity = matrix_equal(m, MATRIX_IDENTITY);
    }
}

void
default_pattern_constructor(Pattern res)
{
    res->ref_count = 0;

    pattern_set_transform(res, MATRIX_IDENTITY);

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
    color_copy(res->fields.concrete.a, a);
    color_copy(res->fields.concrete.b, b);

    res->pattern_at = checker_pattern_at;
}

void
gradient_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = GRADIENT_PATTERN;
    color_copy(res->fields.concrete.a, a);
    color_copy(res->fields.concrete.b, b);

    res->pattern_at = gradient_pattern_at;
}

void
radial_gradient_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = RADIAL_GRADIENT_PATTERN;
    color_copy(res->fields.concrete.a, a);
    color_copy(res->fields.concrete.b, b);

    res->pattern_at = radial_gradient_pattern_at;
}

void
ring_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = RING_PATTERN;
    color_copy(res->fields.concrete.a, a);
    color_copy(res->fields.concrete.b, b);

    res->pattern_at = ring_pattern_at;
}

void
stripe_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);

    res->type = STRIPE_PATTERN;
    color_copy(res->fields.concrete.a, a);
    color_copy(res->fields.concrete.b, b);

    res->pattern_at = stripe_pattern_at;
}

void
uv_align_check_pattern(Color main, Color ul, Color ur, Color bl, Color br, Pattern res)
{
    default_pattern_constructor(res);
    res->type = UV_ALIGN_CHECKER_PATTERN;

    color_copy(res->fields.uv_align_check.main, main);
    color_copy(res->fields.uv_align_check.ul, ul);
    color_copy(res->fields.uv_align_check.ur, ur);
    color_copy(res->fields.uv_align_check.bl, bl);
    color_copy(res->fields.uv_align_check.br, br);

    res->uv_pattern_at = uv_align_check_uv_pattern_at;
}

void
uv_check_pattern(Color a, Color b, size_t width, size_t height, Pattern res)
{
    default_pattern_constructor(res);
    res->type = UV_CHECKER_PATTERN;

    color_copy(res->fields.uv_check.a, a);
    color_copy(res->fields.uv_check.b, b);
    res->fields.uv_check.width = width;
    res->fields.uv_check.height = height;

    res->uv_pattern_at = uv_check_uv_pattern_at;
}

void
uv_gradient_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);
    res->type = UV_GRADIENT_PATTERN;

    color_copy(res->fields.concrete.a, a);
    color_copy(res->fields.concrete.b, b);

    res->uv_pattern_at = uv_gradient_uv_pattern_at;
}

void
uv_radial_gradient_pattern(Color a, Color b, Pattern res)
{
    default_pattern_constructor(res);
    res->type = UV_RADIAL_GRADIENT_PATTERN;

    color_copy(res->fields.concrete.a, a);
    color_copy(res->fields.concrete.b, b);

    res->uv_pattern_at = uv_radial_gradient_uv_pattern_at;
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
pattern_free(Pattern p)
{
    enum uv_map_type type;
    Pattern faces;
    int i;

    if (p != NULL) {
        switch (p->type) {
        case TEXTURE_MAP_PATTERN:
            type = p->fields.uv_map.type;
            faces = p->fields.uv_map.uv_faces;

            switch (type) {
            case CUBE_UV_MAP:
                for (i = 0; i < 6; ++i) {
                    pattern_free(faces + i);
                }
                break;
            case CYLINDER_UV_MAP:
                for (i = 0; i < 3; ++i) {
                    pattern_free(faces + i);
                }
                break;
            case PLANE_UV_MAP:
                for (i = 0; i < 1; ++i) {
                    pattern_free(faces + i);
                }
                break;
            case SPHERE_UV_MAP:
                for (i = 0; i < 1; ++i) {
                    pattern_free(faces + i);
                }
                break;
            case TOROID_UV_MAP:
                for (i = 0; i < 1; ++i) {
                    pattern_free(faces + i);
                }
                break;
            case TRIANGLE_UV_MAP:
                for (i = 0; i < 1; ++i) {
                    pattern_free(faces + i);
                }
                break;
            default:
                break;
            }
            break;
        case NESTED_PATTERN:
            pattern_free(p->fields.nested.pattern1);
            pattern_free(p->fields.nested.pattern2);
            pattern_free(p->fields.nested.pattern3);
            break;
        case BLENDED_PATTERN:
            pattern_free(p->fields.blended.pattern1);
            pattern_free(p->fields.blended.pattern2);
            break;
        case PERTURBED_PATTERN:
            pattern_free(p->fields.perturbed.pattern1);
            break;
        default:
            break;
        }

        p->ref_count--;
        if (p->ref_count == 0) {
            free(p);
        }
    }
}


void
texture_map_pattern(Pattern faces, enum uv_map_type type, Pattern res)
{
    default_pattern_constructor(res);
    res->type = TEXTURE_MAP_PATTERN;

    res->fields.uv_map.type = type;
    res->fields.uv_map.uv_faces = faces;
    res->pattern_at = texture_map_pattern_at;

    int i;

    switch (type) {
    case CUBE_UV_MAP:
        res->uv_map = cube_uv_map;
        for (i = 0; i < 6; ++i) {
            (faces + i)->ref_count += 1;
        }
        break;
    case CYLINDER_UV_MAP:
        res->uv_map = cylinder_uv_map;
        for (i = 0; i < 3; ++i) {
            (faces + i)->ref_count += 1;
        }
        break;
    case PLANE_UV_MAP:
        res->uv_map = plane_uv_map;
        for (i = 0; i < 1; ++i) {
            (faces + i)->ref_count += 1;
        }
        break;
    case SPHERE_UV_MAP:
        res->uv_map = sphere_uv_map;
        for (i = 0; i < 1; ++i) {
            (faces + i)->ref_count += 1;
        }
        break;
    case TOROID_UV_MAP:
        res->uv_map = toroid_uv_map;
        for (i = 0; i < 1; ++i) {
            (faces + i)->ref_count += 1;
        }
        break;
    case TRIANGLE_UV_MAP:
        res->uv_map = triangle_uv_map;
        for (i = 0; i < 1; ++i) {
            (faces + i)->ref_count += 1;
        }
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
texture_map_pattern_alloc(Pattern faces,  enum uv_map_type type)
{
    Pattern p = (Pattern) malloc(sizeof(struct pattern));
    texture_map_pattern(faces, type, p);
    return p;
}
