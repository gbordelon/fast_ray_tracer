#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "renderer.h"
#include "thpool.h"

#include "linalg.h"
#include "shapes.h"
#include "sphere.h"
#include "plane.h"
#include "cube.h"
#include "cone.h"
#include "cylinder.h"
#include "triangle.h"
#include "csg.h"
#include "group.h"
#include "toroid.h"

void
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, Color res);

double
jitter_by(bool jitter)
{
    if (jitter) {
        return (double)rand() / (double)RAND_MAX;
    }
    return 0.5;
}

void
area_light_point_on_light(Light l, size_t u, size_t v, double retval[4])
{
    double uvec[4] = {
        l->u.area.uvec[0],
        l->u.area.uvec[1],
        l->u.area.uvec[2],
        l->u.area.uvec[3]
    };

    double vvec[4] = {
        l->u.area.vvec[0],
        l->u.area.vvec[1],
        l->u.area.vvec[2],
        l->u.area.vvec[3]
    };

    array_scale(uvec, u + jitter_by(l->u.area.jitter));
    array_scale(vvec, v + jitter_by(l->u.area.jitter));

    retval[0] = l->u.area.corner[0] + uvec[0] + vvec[0];
    retval[1] = l->u.area.corner[1] + uvec[1] + vvec[1];
    retval[2] = l->u.area.corner[2] + uvec[2] + vvec[2];
}

#define AREA_LIGHT_CACHE_SIZE 1024
Points
area_light_surface_points(Light light)
{
    if (light->surface_points_cache == NULL) {
        int u, v, i;
        double point[4] = {0.0, 0.0, 0.0, 1.0};
        Points pts = (Points) malloc(AREA_LIGHT_CACHE_SIZE * sizeof(struct pts));
        Points itr = pts;
        for (i = 0, itr = pts; i < AREA_LIGHT_CACHE_SIZE; i++, itr++) {
            itr->points_num = light->num_samples;
            itr->points = (Point) malloc(light->num_samples * sizeof(struct pt));
            for (v = 0; v < light->u.area.vsteps; v++) {
                for (u = 0; u < light->u.area.usteps; u++) {
                    area_light_point_on_light(light, u, v, point);
                    memcpy((itr->points + v * light->u.area.usteps + u)->arr, point, 4 * sizeof(double));
                }
            }
        }
        light->surface_points_cache = pts;
        light->surface_points_cache_len = AREA_LIGHT_CACHE_SIZE;
    }

    int choice = rand() % light->surface_points_cache_len;
    return light->surface_points_cache + choice;
}

Points
point_light_surface_points(Light light)
{
    if (light->surface_points_cache == NULL) {
        Points pts = (Points) malloc(sizeof(struct pts));
        pts->points_num = 1;
        pts->points = (Point) malloc(pts->points_num * sizeof(struct pt));
        memcpy(&pts->points[0].arr, light->u.point.position, sizeof(pts->points[0].arr));
        light->surface_points_cache = pts;
        light->surface_points_cache_len = 1;
    }
    return light->surface_points_cache;
}

double
area_light_intensity_at(Light light, World w, Point p)
{
    double total = 0.0;
    int i;
    Points pts = light->light_surface_points(light);
    Point pos;
    for (i = 0, pos = pts->points; i < pts->points_num; i++, pos++) {
        if (!is_shadowed(w, pos->arr, p)) {
            total += 1.0;
        }
    }

    return total / light->num_samples;
}

double
point_light_intensity_at(Light light, World w, Point p)
{
    if (is_shadowed(w, light->u.point.position, p)) {
        return 0.0;
    }
    return 1.0;
}

void
area_light(double corner[4]/*point*/,
           double full_uvec[4]/*vector*/,
           size_t usteps,
           double full_vvec[4]/*vector*/,
           size_t vsteps,
           bool jitter,
           double intensity[4],
           Light l)
{
    l->type = AREA_LIGHT;

    memcpy(l->u.area.corner, corner, 4 * sizeof(double));

    memcpy(l->u.area.uvec, full_uvec, 4 * sizeof(double));
    array_scale(l->u.area.uvec, 1.0 / (double) usteps);
    l->u.area.usteps = usteps;

    memcpy(l->u.area.vvec, full_vvec, 4 * sizeof(double));
    array_scale(l->u.area.vvec, 1.0 / (double) vsteps);
    l->u.area.vsteps = vsteps;
    l->u.area.jitter = jitter;

    memcpy(l->intensity, intensity, sizeof(l->intensity));
    l->num_samples = usteps * vsteps;


    l->light_surface_points = area_light_surface_points;
    l->intensity_at = area_light_intensity_at;

    // populate surface_points_cache
    l->surface_points_cache = NULL;
    area_light_surface_points(l);
}

Light
area_light_alloc(double corner[4]/*point*/,
                 double full_uvec[4]/*vector*/,
                 size_t usteps,
                 double full_vvec[4]/*vector*/,
                 size_t vsteps,
                 bool jitter,
                 double intensity[4])
{
    Light l = (Light) malloc(sizeof(struct light));
    area_light(corner, full_uvec, usteps, full_vvec, vsteps, jitter, intensity, l);
    return l;
}


void
point_light(Point p, Color intensity, Light l)
{
    l->type = POINT_LIGHT;
    memcpy(l->intensity, intensity, 3 * sizeof(double));
    l->num_samples = 1;
    memcpy(l->u.point.position, p->arr, 4 * sizeof(double));
    l->light_surface_points = point_light_surface_points;
    l->intensity_at = point_light_intensity_at;


    // populate surface_points_cache
    l->surface_points_cache = NULL;
    point_light_surface_points(l);
}

Light
point_light_alloc(Point p, Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    // null check l
    point_light(p, intensity, l);
    return l;
}

Light
array_of_lights(size_t num)
{
    return (Light)malloc(num * sizeof(struct light));
}

Camera
camera(size_t hsize,
       size_t vsize,
       double field_of_view,
       double canvas_distance,
       struct aperture aperture,
       size_t sample_num,
       Matrix transform)
{
    Camera c = (Camera) malloc(sizeof(struct camera));
    // null check c

    c->hsize = hsize;
    c->vsize = vsize;
    c->field_of_view = field_of_view;
    c->canvas_distance = canvas_distance;
    c->sample_num = sample_num;
    c->aperture = aperture;


    double half_view = canvas_distance * tan(field_of_view * 0.5);
    double aspect = (double)hsize / (double)vsize;

    if (aspect >= 1.0) {
        c->half_width = half_view;
        c->half_height = half_view / aspect;
    } else {
        c->half_width = half_view * aspect;
        c->half_height = half_view;
    }

    c->transform = NULL;
    c->transform_inverse = NULL;
    camera_set_transform(c, transform);
    c->pixel_size = c->half_width * 2.0 / (double)hsize;

    return c;
}

void
camera_set_transform(Camera c, Matrix m)
{
    if (c != NULL) {
        if (c->transform_inverse != NULL) {
            matrix_free(c->transform_inverse);
        }
        if (c->transform != NULL) {
            matrix_free(c->transform);
        }
        c->transform = m;
        c->transform_inverse = matrix_inverse_alloc(m);
    }
}

Matrix
view_transform(Point fr, Point to, Vector up)
{
    struct v v;
    struct v forward;
    struct v upn;
    struct v left;
    struct v true_up;
    struct m orientation;
    struct m m;

    vector_from_points(to,fr, &v);
    vector_normalize(&v, &forward);

    vector_normalize(up, &upn);
    vector_cross(&forward, &upn, &left);
    vector_cross(&left, &forward, &true_up);

    matrix(left.arr[0], left.arr[1], left.arr[2], 0,
           true_up.arr[0], true_up.arr[1], true_up.arr[2], 0,
           -forward.arr[0], -forward.arr[1], -forward.arr[2], 0,
           0, 0, 0, 1,
           &orientation);

    matrix_translate(-fr->arr[0], -fr->arr[1], -fr->arr[2], &m);

    return matrix_multiply_alloc(&orientation, &m);
}

World
world()
{
    World w = (World) malloc(sizeof(struct world));
    w->lights = NULL;
    w->lights_num = 0;
    w->shapes = NULL;
    w->shapes_num = 0;
    w->xs = intersections_empty(64);
    return w;
}

World
default_world()
{
    World w = world();
    Point p = point(-10, 10, -10);
    Color c = color(1.0, 1.0, 1.0);
    Light l = point_light_alloc(p, c);
    w->lights = l;
    w->lights_num = 1;
    Shape shapes = (Shape) malloc(2 * sizeof(struct shape));

    Shape s1 = shapes;
    Shape s2 = shapes + 1;

    toroid(s1);
    sphere(s2);

    s1->material->color[0] = 0.8;
    s1->material->color[1] = 1.0;
    s1->material->color[2] = 0.6;
    s1->material->diffuse = 0.7;
    s1->material->specular = 0.2;
    s1->material->casts_shadow = true;
    s1->material->transparency = 0.0;
    s1->fields.toroid.r1 = 1.0;
    s1->fields.toroid.r2 = 0.5;

    Matrix scaling = matrix_scale_alloc(0.5, 0.5, 0.5);
    shape_set_transform(s2, scaling);
    shape_set_transform(s1, transform_chain(matrix_rotate_x_alloc(M_PI_4), matrix_scale_alloc(3,3,5)));

    w->shapes = shapes;
    w->shapes_num = 2;

    return w;
}

bool
is_shadowed(World w, double light_position[4], Point pt)
{
    struct v v;
    vector_from_arrays(light_position, pt->arr, &v); // from pt to light

    double distance = vector_magnitude(&v); // distance between pt and light

    struct v direction;
    vector_normalize(&v, &direction);

    struct ray r;
    ray_array(pt->arr, direction.arr, &r);

    Intersections xs = intersect_world(w, &r);

    Intersection h = hit(xs, true);
    bool retval = h != NULL && h->t < distance;

    return retval;
}


void
circle_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double u, v;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
    } while (u * u + v * v > bounds->u.circle.r1);
}

void
cross_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double u, v;
    bool check;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        check = ((u > bounds->u.cross.x1) && (u <= bounds->u.cross.x2))
             || ((v > bounds->u.cross.y1) && (v <= bounds->u.cross.y2));
    } while (!check);
}

void
diamond_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double u, v;
    bool check;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        check = (u <= 0)
              ? (-u + bounds->u.diamond.b1 <= v) && (v < u + bounds->u.diamond.b2)
              : (0 <= *x)
              ? (u + bounds->u.diamond.b3 <= v) && (v < -u + bounds->u.diamond.b4)
              : false;
    } while (!check);
}

void
double_circle_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double mag;
    double u, v;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        mag = u * u + v * v;
    } while (mag > bounds->u.double_circle.r1 || mag < bounds->u.double_circle.r2);
}

void
point_aperture_fn(double *x, double *y, struct aperture *bound)
{
    *x = 0.5;
    *x = 0.5;
}

void
square_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    *x = jitter_by(true);
    *y = jitter_by(true);
}

void
random_point_by_function(double *x, double *y, struct aperture *bounds, void (*fn)(double *, double *, struct aperture *))
{
    fn(x, y, bounds);
    *x -= 0.5;
    *x -= 0.5;
}

void
ray_for_pixel(Camera cam, double px, double py, double x_jitter, double y_jitter, Ray res)
{
    double xoffset = (px + x_jitter) * cam->pixel_size;
    double yoffset = (py + y_jitter) * cam->pixel_size;
    double world_x = cam->half_width - xoffset;
    double world_y = cam->half_height - yoffset;
    Matrix inv = cam->transform_inverse;
    struct pt origin;
    struct v direction;
    struct pt p;
    struct pt pixel;
    struct v v;
    void (*aperture_fn)(double *, double *, struct aperture *);

    switch (cam->aperture.type) {
    case CIRCULAR_APERTURE:
        aperture_fn = circle_aperture_fn;
        break;
    case CROSS_APERTURE:
        aperture_fn = cross_aperture_fn;
        break;
    case DIAMOND_APERTURE:
        aperture_fn = diamond_aperture_fn;
        break;
    case DOUBLE_CIRCLE_APERTURE:
        aperture_fn = double_circle_aperture_fn;
        break;
    case SQUARE_APERTURE:
        aperture_fn = square_aperture_fn;
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
        aperture_fn = point_aperture_fn;
        break;
    }

    p.arr[0] = world_x;
    p.arr[1] = world_y;
    p.arr[2] = -cam->canvas_distance;
    p.arr[3] = 1.0;
    matrix_point_multiply(inv, &p, &pixel);

    p.arr[0] = 0;
    p.arr[1] = 0;
    random_point_by_function(p.arr, p.arr+1, &cam->aperture, aperture_fn);

    p.arr[0] *= cam->aperture.size;
    p.arr[1] *= cam->aperture.size;
    p.arr[2] = 0;
    matrix_point_multiply(inv, &p, &origin);
    vector_from_points(&pixel, &origin, &v);
    vector_normalize(&v, &direction);

    ray_array(origin.arr, direction.arr, res);
}

/*
 * ^  +-----+
 * |  |     |
 * |  |     |
 * |  +-----+
 * v u------>
 *
 * u ranges from x to x+1
 * v ranges from y to y+1
 * subdivide the pixel by units to make a grid
 * choose a point in each grid cell to sample
 * average the colors
 */
void
pixel_multi_sample(Camera cam, World w, double x, double y, size_t usteps, size_t vsteps, bool jitter, Color res)
{
    double x_offset, y_offset;
    size_t u, v;
    double total_steps = (double)usteps * (double)vsteps;
    double usteps_inv = 1.0 / (double)usteps;
    double vsteps_inv = 1.0 / (double)vsteps;
    struct ray r;
    struct color c;
    struct color acc;

    acc.arr[0] = 0;
    acc.arr[1] = 0;
    acc.arr[2] = 0;

    for (v = 0; v < vsteps; v++) {
        for (u = 0; u < usteps; u++) {
            x_offset = ((double)u + jitter_by(jitter)) * usteps_inv;
            y_offset = ((double)v + jitter_by(jitter)) * vsteps_inv;
            c.arr[0] = 0;
            c.arr[1] = 0;
            c.arr[2] = 0;
            r.origin[0] = 0;
            r.origin[1] = 0;
            r.origin[2] = 0;
            r.origin[3] = 1;
            r.direction[0] = 0;
            r.direction[1] = 0;
            r.direction[2] = 0;
            r.direction[0] = 0;
            ray_for_pixel(cam, x, y, x_offset, y_offset, &r);
            color_at(w, &r, 5, &c);
            color_accumulate(&acc, &c);
        }
    }

    color_scale(&acc, 1.0 / total_steps);
    memcpy(res->arr, acc.arr, 3 * sizeof(double));
}

Canvas
render(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter)
{
    int i,j,k;
    struct color c;

    Canvas image = canvas_alloc(cam->hsize, cam->vsize);

    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            c.arr[0] = 0;
            c.arr[1] = 0;
            c.arr[2] = 0;
            pixel_multi_sample(cam, w, (double)i, (double)j, usteps, vsteps, jitter, &c);
            canvas_write_pixel(image, i, j, &c);
        }
        k += 1;
        printf("Wrote %d rows out of %lu\n", k, cam->vsize);
    }

    return image;
}

struct render_args {
    Camera cam;
    World w;
    Canvas image;
    size_t x;
    size_t y;
    size_t usteps;
    size_t vsteps;
    bool jitter;
};

void
render_multi_helper(void *args)
{
    Camera cam = ((struct render_args *)args)->cam;
    World w = ((struct render_args *)args)->w;
    Canvas image = ((struct render_args *)args)->image;
    size_t x = ((struct render_args *)args)->x;
    size_t y = ((struct render_args *)args)->y;
    size_t usteps = ((struct render_args *)args)->usteps;
    size_t vsteps = ((struct render_args *)args)->vsteps;
    bool jitter = ((struct render_args *)args)->jitter;
    struct color c;

    c.arr[0] = 0;
    c.arr[1] = 0;
    c.arr[2] = 0;

    pixel_multi_sample(cam, w, (double)x, (double)y, usteps, vsteps, jitter, &c);
    canvas_write_pixel(image, x, y, &c);
}

/*
 * World is not thread safe because of Intersections buffers
 *
 * Use realloc instead of malloc to keep references in tact when a thread needs to
 * resize an Intersections array
 *
 * 
 */
Canvas
render_multi(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter)
{
    int i,j,k;
    threadpool thpool = thpool_init(1);
    Canvas image = canvas_alloc(cam->hsize, cam->vsize);

    struct render_args *args_array = (struct render_args *)malloc(cam->vsize * cam->hsize * sizeof(struct render_args));
    
    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            (args_array + k)->cam = cam;
            (args_array + k)->w = w;
            (args_array + k)->image = image;
            (args_array + k)->x = i;
            (args_array + k)->y = j;
            (args_array + k)->usteps = usteps;
            (args_array + k)->vsteps = vsteps;
            (args_array + k)->jitter = jitter;
            thpool_add_work(thpool, render_multi_helper, (void*)(args_array + k));
            k++;
        }
        printf("adding work to thpool %d\n", k);
    }

    thpool_wait(thpool);

    thpool_destroy(thpool);
    free(args_array);

    return image;
}

void
color_at(World w, Ray r, size_t remaining, Color res)
{
    Intersections xs = intersect_world(w, r);
    Intersection i = hit(xs, false);
    struct computations comps;
    struct color c;

    c.arr[0] = 0;
    c.arr[1] = 0;
    c.arr[2] = 0;

    if (i != NULL) {
        prepare_computations(i, r, xs, &comps);
        shade_hit(w, &comps, remaining, &c);
    }
    memcpy(res->arr, c.arr, 3 * sizeof(double));
}

int
sort_intersections_asc(const void *p, const void *q)
{
    double l = ((Intersection)p)->t;
    double r = ((Intersection)q)->t;
    if (l - r < 0) {
        return -1;
    } else if (l - r > 0) {
        return 1;
    }
    return 0;
}

int
sort_intersections_desc(const void *p, const void *q)
{
    double l = ((Intersection)p)->t;
    double r = ((Intersection)q)->t;
    if (l - r < 0) {
        return 1;
    } else if (l - r > 0) {
        return -1;
    }
    return 0;
}

void
intersections_reverse(Intersections xs)
{
    // sort xs by xs->xs->t descending
    mergesort((void*)xs->xs, xs->num, sizeof(struct intersection), sort_intersections_desc);
}

void
intersections_sort(Intersections xs)
{
    // sort xs by xs->xs->t ascending
    mergesort((void*)xs->xs, xs->num, sizeof(struct intersection), sort_intersections_asc);
}

Intersections
intersect_world(World w, Ray r)
{
    int i;

    w->xs->num = 0;
    for (i = 0; i < w->shapes_num; i++) {
        Intersections xs_1 = (w->shapes + i)->intersect(w->shapes + i, r);
        if (xs_1 == NULL || xs_1->num == 0) {
            continue;
        }

        // realloc
        if (xs_1->num + w->xs->num >= w->xs->array_len) {
            Intersections xs_2 = intersections_empty(2 * w->xs->array_len);
            memcpy(xs_2->xs, w->xs->xs, w->xs->array_len * sizeof(struct intersection));
            xs_2->num = w->xs->num;
            intersections_free(w->xs);
            w->xs = xs_2;
        }

        // copy from xs_1 into xs + xs->num
        memcpy(w->xs->xs + w->xs->num, xs_1->xs, xs_1->num * sizeof(struct intersection));
        w->xs->num += xs_1->num;
    }

    if (w->xs->num > 1) {
        intersections_sort(w->xs);
    }

    return w->xs;
}

void
position(Ray ray, double t, Point position)
{
    memcpy(position->arr, ray->origin, 4 * sizeof(double));
    position->arr[0] += ray->direction[0] * t;
    position->arr[1] += ray->direction[1] * t;
    position->arr[2] += ray->direction[2] * t;
}

void
prepare_computations(Intersection i, Ray r, Intersections xs, Computations res)
{
    res->t = i->t;

    res->obj = i->object;

    position(r, i->t, &res->p);

    i->object->normal_at(i->object, &res->p, i, &res->normalv);

    memcpy(&res->eyev.arr, r->direction, 4 * sizeof(double));
    array_scale(res->eyev.arr, -1.0);

    res->inside = false;

    vector_array_reflect(r->direction, res->normalv.arr, &res->reflectv);

    if (array_dot(res->normalv.arr, res->eyev.arr) < 0) {
        res->inside = true;
        array_scale(res->normalv.arr, -1);
    }

    res->over_point.arr[0] = res->p.arr[0] + res->normalv.arr[0] * EPSILON;
    res->over_point.arr[1] = res->p.arr[1] + res->normalv.arr[1] * EPSILON;
    res->over_point.arr[2] = res->p.arr[2] + res->normalv.arr[2] * EPSILON;
    res->over_point.arr[3] = 1.0;

    res->under_point.arr[0] = res->p.arr[0] - res->normalv.arr[0] * EPSILON;
    res->under_point.arr[1] = res->p.arr[1] - res->normalv.arr[1] * EPSILON;
    res->under_point.arr[2] = res->p.arr[2] - res->normalv.arr[2] * EPSILON;
    res->under_point.arr[3] = 1.0;

    res->n1 = 1.0;
    res->n2 = 1.0;

    Computations c = res;

    int j, k;
    Intersection x;
    size_t container_len = 0;
    Shape *container = (Shape*) malloc(xs->num * sizeof(Shape));

    for (j = 0, x = xs->xs; j < xs->num; x++, j++) {
        if (x == i) {// address compare should be okay.
            if (container_len == 0) {
                c->n1 = 1.0;
            } else {
                c->n1 = container[container_len-1]->material->refractive_index;
            }
        }

        for (k = 0; k < container_len; k++) {
            if (container[k] == i->object) {
                break;
            }
        }

        if (k < container_len) {
            // shift everything left one slot, overwriting container[index_of_object] first
            for (; k < container_len - 1; k++) {
                container[k] = container[k+1];
            }
            container_len--;
        } else {
            container[container_len] = i->object;
            container_len++;
        }

        if (x == i) {
            if (container_len == 0) {
                c->n2 = 1.0;
            } else {
                c->n2 = container[container_len-1]->material->refractive_index;
            }
            break;
        }
    }

    free(container);
}

void
reflected_color(World w, Computations comps, size_t remaining, Color res)
{
    if (remaining == 0 || equal(comps->obj->material->reflective, 0)) {
        res->arr[0] = 0;
        res->arr[1] = 0;
        res->arr[2] = 0;
    } else {
        struct ray reflect_ray;
        struct color c;
        c.arr[0] = 0;
        c.arr[1] = 0;
        c.arr[2] = 0;
        ray_array(comps->over_point.arr, comps->reflectv.arr, &reflect_ray);
        color_at(w, &reflect_ray, remaining - 1, &c);
        color_scale(&c, comps->obj->material->reflective);
        color_accumulate(res, &c);
    }
}

void
refracted_color(World w, Computations comps, size_t remaining, Color res)
{
    if (remaining == 0 || equal(comps->obj->material->transparency,0)) {
        res->arr[0] = 0;
        res->arr[1] = 0;
        res->arr[2] = 0;
        return;
    }

    double n_ratio = comps->n1 / comps->n2;
    double cos_i = vector_dot(&comps->eyev, &comps->normalv);
    double sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);

    if (sin2_t > 1.0) {
        res->arr[0] = 0;
        res->arr[1] = 0;
        res->arr[2] = 0;
        return;
    }

    struct color c;
    c.arr[0] = 0;
    c.arr[1] = 0;
    c.arr[2] = 0;

    struct v t1;
    struct v t2;
    struct v direction;
    struct ray refracted_ray;

    double cos_t = sqrt(1.0 - sin2_t);
    memcpy(t1.arr, comps->normalv.arr, 4 * sizeof(double));
    vector_scale(&t1, n_ratio * cos_i - cos_t);
    memcpy(t2.arr, comps->eyev.arr, 4 * sizeof(double));

    vector_scale(&t2, n_ratio);
    direction.arr[0] = t1.arr[0] - t2.arr[0];
    direction.arr[1] = t1.arr[1] - t2.arr[1];
    direction.arr[2] = t1.arr[2] - t2.arr[2];
    direction.arr[3] = 0.0;

    ray_array(comps->under_point.arr, direction.arr, &refracted_ray);
    color_at(w, &refracted_ray, remaining - 1, &c);
    color_scale(&c, comps->obj->material->transparency);
    color_accumulate(res, &c);
}

double
schlick(Computations comps)
{
    double co = vector_dot(&comps->eyev, &comps->normalv);
    if (comps->n1 > comps->n2) {
        double n = comps->n1 / comps->n2;
        double sin2_t = n * n * (1.0 - co * co);
        if (sin2_t > 1.0) {
            return 1.0;
        }
        co = sqrt(1.0 - sin2_t);
    }

    double r0 = ((comps->n1 - comps->n2) / (comps->n1 + comps->n2));
    r0 = r0 * r0;

    return r0 + (1.0 - r0) * (1.0 - co) * (1.0 - co) * (1.0 - co) * (1.0 - co) * (1.0 - co);
}

void
shade_hit(World w, Computations comps, size_t remaining, Color res)
{
    Light itr;
    size_t i;
    struct color c;
    struct color surface;
    double intensity;

    surface.arr[0] = 0;
    surface.arr[1] = 0;
    surface.arr[2] = 0;

    for (i = 0, itr = w->lights; i < w->lights_num; i++, itr++) {
        c.arr[0] = 0;
        c.arr[1] = 0;
        c.arr[2] = 0;
        intensity = itr->intensity_at(itr, w, &comps->over_point);
        lighting(comps->obj->material,
                 comps->obj,
                 itr,
                 &comps->over_point,
                 &comps->eyev,
                 &comps->normalv,
                 intensity,
                 &c);

        color_accumulate(&surface, &c);
    }

    struct color reflected;
    reflected.arr[0] = 0;
    reflected.arr[1] = 0;
    reflected.arr[2] = 0;

    reflected_color(w, comps, remaining, &reflected);

    struct color refracted;
    refracted.arr[0] = 0;
    refracted.arr[1] = 0;
    refracted.arr[2] = 0;

    refracted_color(w, comps, remaining, &refracted);

    if (comps->obj->material->reflective > 0 && comps->obj->material->transparency > 0) {
        double reflectance = schlick(comps);
        color_scale(&reflected, reflectance);
        color_scale(&refracted, 1.0 - reflectance);
    }

    color_accumulate(&surface, &reflected);
    color_accumulate(&surface, &refracted);

    memcpy(res->arr, surface.arr, 3 * sizeof(double));
}

void
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, Color res)
{
    struct color scolor;
    struct color ambient;

    if (material->pattern != NULL) {
        material->pattern->pattern_at_shape(material->pattern, shape, point, &scolor);
    } else {
        memcpy(scolor.arr, material->color, 3 * sizeof(double));
    }

    scolor.arr[0] *= light->intensity[0];
    scolor.arr[1] *= light->intensity[1];
    scolor.arr[2] *= light->intensity[2];

    memcpy(ambient.arr, scolor.arr, 3 * sizeof(double));

    color_scale(&ambient, material->ambient);
    if (equal(shade_intensity, 0.0)) {
        color_accumulate(res, &ambient);
        return;
    }

    struct color diffuse;
    int i;
    Points pts = light->light_surface_points(light);
    Point position;
    struct v diff;
    struct v lightv;
    struct v reflectv;

    for (i = 0, position = pts->points; i < pts->points_num; i++, position++) {
        vector_from_points(position, point, &diff);
        vector_normalize(&diff, &lightv);
        double light_dot_normal = vector_dot(&lightv, normalv);
        if (light_dot_normal >= 0) {
            // diffuse
            memcpy(diffuse.arr, scolor.arr, 3 * sizeof(double));
            color_scale(&diffuse, material->diffuse);
            color_scale(&diffuse, light_dot_normal);
            color_accumulate(res, &diffuse);

            // specular
            vector_scale(&lightv, -1);
            vector_reflect(&lightv, normalv, &reflectv);
            vector_scale(&lightv, -1);

            double reflect_dot_eye = vector_dot(&reflectv, eyev);
            if (reflect_dot_eye > 0 && material->specular > 0) {
                double factor = pow(reflect_dot_eye, material->shininess);
                res->arr[0] += light->intensity[0] * material->specular * factor;
                res->arr[1] += light->intensity[1] * material->specular * factor;
                res->arr[2] += light->intensity[2] * material->specular * factor;
            }
        }
    }

    double scaling = shade_intensity / (double)light->num_samples;
    color_scale(res, scaling);
    color_accumulate(res, &ambient);
}
