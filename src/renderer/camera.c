#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libs/linalg/linalg.h"
#include "../libs/sampler/sampler.h"

#include "camera.h"


void
circle_aperture_fn(double *x, double *y, Aperture aperture)
{
    double u, v;
    do {
        *x = drand48();
        *y = drand48();
        u = 2 * *x - 1;
        v = 2 * *y - 1;
    } while (u * u + v * v > aperture->u.circle.r1);
}

void
cross_aperture_fn(double *x, double *y, Aperture aperture)
{
    double u, v;
    bool check;
    do {
        *x = drand48();
        *y = drand48();
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        check = ((u > aperture->u.cross.x1) && (u <= aperture->u.cross.x2))
             || ((v > aperture->u.cross.y1) && (v <= aperture->u.cross.y2));
    } while (!check);
}

void
diamond_aperture_fn(double *x, double *y, Aperture aperture)
{
    double u, v;
    bool check;
    do {
        *x = drand48();
        *y = drand48();
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        check = (u <= 0)
              ? (-u + aperture->u.diamond.b1 <= v) && (v < u + aperture->u.diamond.b2)
              : (0 <= *x)
              ? (u + aperture->u.diamond.b3 <= v) && (v < -u + aperture->u.diamond.b4)
              : false;
    } while (!check);
}

void
doughnut_aperture_fn(double *x, double *y, Aperture aperture)
{
    double mag;
    double u, v;
    do {
        *x = drand48();
        *y = drand48();
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        mag = u * u + v * v;
    } while (mag > aperture->u.doughnut.r1 || mag < aperture->u.doughnut.r2);
}

void
point_aperture_fn(double *x, double *y, Aperture aperture)
{
    *x = 0.5;
    *y = 0.5;
}

void
square_aperture_fn(double *x, double *y, Aperture aperture)
{
    *x = drand48();
    *y = drand48();
}

void
sample_aperture(double xy[2], size_t u, size_t v, const Aperture aperture)
{
    aperture->aperture_fn(xy, xy+1, aperture);
    *xy -= 0.5;
    *(xy+1) -= 0.5;
}


void
camera_set_transform(Camera c, Matrix m)
{
    if (c != NULL) {
        matrix_copy(m, c->transform);
        matrix_inverse(m, c->transform_inverse);
    }
}

Camera
camera(size_t hsize,
       size_t vsize,
       double field_of_view,
       double canvas_distance,
       size_t usteps,
       size_t vsteps,
       Aperture aperture,
       Matrix transform)
{
    Camera c = (Camera) malloc(sizeof(struct camera));

    c->hsize = hsize;
    c->vsize = vsize;
    c->field_of_view = field_of_view;
    c->canvas_distance = canvas_distance;
    c->aperture = *aperture;
    c->usteps = usteps;
    c->vsteps = vsteps;

    camera_set_transform(c, transform);

    double half_view = canvas_distance * tan(field_of_view * 0.5);
    double aspect = (double)hsize / (double)vsize;

    if (aspect >= 1.0) {
        c->half_width = half_view;
        c->half_height = half_view / aspect;
    } else {
        c->half_width = half_view * aspect;
        c->half_height = half_view;
    }

    c->pixel_size = c->half_width * 2.0 / (double)hsize;

    return c;
}

void
view_transform(Point fr, Point to, Vector up, Matrix res)
{
    Vector v;
    Vector forward;
    Vector upn;
    Vector left;
    Vector true_up;
    Matrix orientation;
    Matrix m;

    vector_from_points(to, fr, v);
    vector_normalize(v, forward);

    vector_normalize(up, upn);
    vector_cross(forward, upn, left);
    vector_cross(left, forward, true_up);

    matrix(left[0], left[1], left[2], 0,
           true_up[0], true_up[1], true_up[2], 0,
           -forward[0], -forward[1], -forward[2], 0,
           0, 0, 0, 1,
           orientation);

    matrix_translate(-fr[0], -fr[1], -fr[2], m);

    matrix_multiply(orientation, m, res);
}

void
aperture(enum aperture_type type, double size, size_t usteps, size_t vsteps, bool jitter, Aperture res)
{
    res->type = type;
    res->size = size;
    res->jitter = jitter;
    sampler_2d(jitter, usteps, vsteps, sampler_default_constraint, &res->sampler);

    switch (type) {
    case CIRCULAR_APERTURE:
        res->aperture_fn = circle_aperture_fn;
        break;
    case CROSS_APERTURE:
        res->aperture_fn = cross_aperture_fn;
        break;
    case DIAMOND_APERTURE:
        res->aperture_fn = diamond_aperture_fn;
        break;
    case DOUGHNUT_APERTURE:
        res->aperture_fn = doughnut_aperture_fn;
        break;
    case SQUARE_APERTURE:
        res->aperture_fn = square_aperture_fn;
        break;
    case HEXAGONAL_APERTURE:
        // not yet impl
    case PENTAGONAL_APERTURE:
        // not yet impl
    case OCTAGONAL_APERTURE:
        // not yet impl
    case POINT_APERTURE:
    default:
        // no jittering
        // no need to do any focal_blur stuff
        res->aperture_fn = point_aperture_fn;
        break;
    }
}

void
circle_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct circle_aperture_args *args, Aperture res)
{
    aperture(CIRCULAR_APERTURE, size, usteps, vsteps, jitter, res);
    memcpy(&res->u.circle, args, sizeof(struct circle_aperture_args));
}

void cross_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct cross_aperture_args *args, Aperture res)
{
    aperture(CROSS_APERTURE, size, usteps, vsteps, jitter, res);
    memcpy(&res->u.cross, args, sizeof(struct cross_aperture_args));
}

void diamond_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct diamond_aperture_args *args, Aperture res)
{
    aperture(DIAMOND_APERTURE, size, usteps, vsteps, jitter, res);
    memcpy(&res->u.diamond, args, sizeof(struct diamond_aperture_args));
}

void doughnut_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct doughnut_aperture_args *args, Aperture res)
{
    aperture(DOUGHNUT_APERTURE, size, usteps, vsteps, jitter, res);
    memcpy(&res->u.doughnut, args, sizeof(struct doughnut_aperture_args));
}

void square_aperture(double size, size_t usteps, size_t vsteps, bool jitter, Aperture res)
{
    aperture(SQUARE_APERTURE, size, usteps, vsteps, jitter, res);
}

void point_aperture(Aperture res)
{
    aperture(POINT_APERTURE, 0, 1, 1, false, res);
}
