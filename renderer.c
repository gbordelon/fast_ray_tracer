#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "renderer.h"

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

Color
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity);

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
    memcpy(l->intensity, intensity, sizeof(l->intensity));
    l->num_samples = 1;
    memcpy(l->u.point.position, p->arr, sizeof(l->u.point.position));
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
       double aperture_size,
       double canvas_distance,
       enum aperture_shape aperture_shape,
       size_t sample_num,
       bool jitter,
       Matrix transform)
{
    Camera c = (Camera) malloc(sizeof(struct camera));
    // null check c

    c->hsize = hsize;
    c->vsize = vsize;
    c->field_of_view = field_of_view;
    c->aperture_size = aperture_size;
    c->canvas_distance = canvas_distance;
    c->aperture_shape = aperture_shape;
    c->sample_num = sample_num;
    c->jitter = jitter;


    double half_view = canvas_distance * tan(field_of_view / 2.0);
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
    Vector v = vector_from_points_alloc(to,fr);
    Vector forward = vector_normalize_alloc(v);

    Vector upn = vector_normalize_alloc(up);
    Vector left = vector_cross_alloc(forward, upn);
    Vector true_up = vector_cross_alloc(left, forward);
    Matrix orientation = matrix_alloc(left->arr[0], left->arr[1], left->arr[2], 0,
                                      true_up->arr[0], true_up->arr[1], true_up->arr[2], 0,
                                      -forward->arr[0], -forward->arr[1], -forward->arr[2], 0,
                                      0, 0, 0, 1);
    Matrix m = matrix_translate_alloc(-fr->arr[0], -fr->arr[1], -fr->arr[2]);

    Matrix retval = matrix_multiply_alloc(orientation, m);

    matrix_free(m);
    matrix_free(orientation);
    vector_free(true_up);
    vector_free(left);
    vector_free(upn);
    vector_free(forward);
    vector_free(v);

    return retval;
}

World
world()
{
    World w = (World) malloc(sizeof(struct world));
    w->lights = NULL;
    w->lights_num = 0;
    w->shapes = NULL;
    w->shapes_num = 0;
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

    sphere(s1);
    sphere(s2);

    s1->material->color[0] = 0.8;
    s1->material->color[1] = 1.0;
    s1->material->color[2] = 0.6;
    s1->material->diffuse = 0.7;
    s1->material->specular = 0.2;

    Matrix scaling = matrix_scale_alloc(0.5, 0.5, 0.5);
    shape_set_transform(s2, scaling);

    w->shapes = shapes;
    w->shapes_num = 2;

    return w;
}

bool
is_shadowed(World w, double light_position[4], Point pt)
{
    Vector v = vector_from_arrays_alloc(light_position, pt->arr); // from pt to light
    double distance = vector_magnitude(v); // distance between pt and light
    Vector direction = vector_normalize_alloc(v);
    Ray r = ray_alloc(pt, direction);
    Intersections xs = intersect_world(w, r);

    Intersection h = hit(xs, true);
    bool retval = h != NULL && h->t < distance;

    intersections_free(xs);
    ray_free(r);
    vector_free(direction);
    vector_free(v);

    return retval;
}

Ray
ray_for_pixel(Camera cam, double px, double py, double x_offset, double y_offset)
{
    double xoffset = (px + x_offset) * cam->pixel_size;
    double yoffset = (py + y_offset) * cam->pixel_size;
    double world_x = cam->half_width - xoffset;
    double world_y = cam->half_height - yoffset;
    Matrix inv = cam->transform_inverse;
    Point p = point(world_x, world_y, -cam->canvas_distance);

    Point pixel = matrix_point_multiply_alloc(inv, p);
    p->arr[0] = 0 + (-0.5 + jitter_by(true)) * cam->aperture_size; // aperture size of 1
    p->arr[1] = 0 + (-0.5 + jitter_by(true)) * cam->aperture_size; // aperture size of 1
    p->arr[2] = 0;
    Point origin = matrix_point_multiply_alloc(inv, p);
    Vector v = vector_from_points_alloc(pixel, origin);
    Vector direction = vector_normalize_alloc(v);

    Ray r = ray_alloc(origin, direction);

    point_free(p);
    vector_free(v);

    return r;
}

/*
Ray
ray_for_pixel_ul(Camera cam, double px, double py, double factor)
{
    return ray_for_pixel(cam, px, py, 0, factor);
}

Ray
ray_for_pixel_ur(Camera cam, double px, double py, double factor)
{
    return ray_for_pixel(cam, px, py, factor, factor);
}

Ray
ray_for_pixel_ll(Camera cam, double px, double py, double factor)
{
    return ray_for_pixel(cam, px, py, 0, 0);
}

Ray
ray_for_pixel_lr(Camera cam, double px, double py, double factor)
{
    return ray_for_pixel(cam, px, py, 0, factor);
}

Ray
ray_for_pixel_m(Camera cam, double px, double py, double factor)
{
    return ray_for_pixel(cam, px, py, 0.5 * factor, 0.5 * factor);
}

Color
pixel_single_sample(Camera cam, World w, int x, int y)
{
    Ray rm = ray_for_pixel_m(cam, (double)x, (double)y, 1.0);
    Color cm = color_at(w, rm, 5);

    ray_free(rm);

    return cm;
}

Color
pixel_multi_sample(Camera cam, World w, double x, double y, double factor)
{
    double color_threshold = 0.9;
    double new_factor = 0.5 * factor;

    Ray rm = ray_for_pixel_m(cam, x, y, factor);
    Color cm = color_at(w, rm, 5);

    Ray rul = ray_for_pixel_ul(cam, x, y, factor);
    Color cul = color_at(w, rul, 5);
    if (fabs(cul->arr[0] - cm->arr[0]) > color_threshold ||
            fabs(cul->arr[1] - cm->arr[1]) > color_threshold ||
            fabs(cul->arr[2] - cm->arr[2]) > color_threshold) {
        color_free(cul);
        cul = pixel_multi_sample(cam, w, x, y + new_factor, new_factor);
    }

    Ray rur = ray_for_pixel_ur(cam, x, y, factor);
    Color cur = color_at(w, rur, 5);
    if (fabs(cur->arr[0] - cm->arr[0]) > color_threshold ||
            fabs(cur->arr[1] - cm->arr[1]) > color_threshold ||
            fabs(cur->arr[2] - cm->arr[2]) > color_threshold) {
        color_free(cur);
        cur = pixel_multi_sample(cam, w, x + new_factor, y + new_factor, new_factor);
    }

    Ray rll = ray_for_pixel_ll(cam, x, y, factor);
    Color cll = color_at(w, rll, 5);
    if (fabs(cll->arr[0] - cm->arr[0]) > color_threshold ||
            fabs(cll->arr[1] - cm->arr[1]) > color_threshold ||
            fabs(cll->arr[2] - cm->arr[2]) > color_threshold) {
        color_free(cll);
        cll = pixel_multi_sample(cam, w, x, y, new_factor);
    }

    Ray rlr = ray_for_pixel_lr(cam, x, y, factor);
    Color clr = color_at(w, rlr, 5);
    if (fabs(clr->arr[0] - cm->arr[0]) > color_threshold ||
            fabs(clr->arr[1] - cm->arr[1]) > color_threshold ||
            fabs(clr->arr[2] - cm->arr[2]) > color_threshold) {
        color_free(clr);
        clr = pixel_multi_sample(cam, w, x + new_factor, y, new_factor);
    }

    cm->arr[0] += cul->arr[0] + cur->arr[0] + cll->arr[0] + clr->arr[0];
    cm->arr[1] += cul->arr[1] + cur->arr[1] + cll->arr[1] + clr->arr[1];
    cm->arr[2] += cul->arr[2] + cur->arr[2] + cll->arr[2] + clr->arr[2];

    color_scale(cm, 0.2); // divide by 5

    color_free(clr);
    ray_free(rlr);
    color_free(cll);
    ray_free(rll);
    color_free(cur);
    ray_free(rur);
    color_free(cul);
    ray_free(rul);
    ray_free(rm);

    return cm;
}
*/


/*
 * ^  +-----+
 * |  |     |
 * |  |     |
 *    +-----+
 * v u  ---->
 *
 * u ranges from x to x+1
 * v ranges from y to y+1
 * subdivide the pixel by units to make a grid
 * choose a point in each grid cell to sample
 * average the colors
 */
Color
pixel_multi_sample(Camera cam, World w, double x, double y, size_t usteps, size_t vsteps, bool jitter)
{
    double sum[3] = {0.0, 0.0, 0.0};
    double x_offset, y_offset;
    size_t u, v;
    double total_steps = (double)usteps * (double)vsteps;

    Ray r;
    Color c;

    for (v = 0; v < vsteps; v++) {
        for (u = 0; u < usteps; u++) {
            x_offset = ((double)u + jitter_by(jitter)) / (double)usteps;
            y_offset = ((double)v + jitter_by(jitter)) / (double)vsteps;
            r = ray_for_pixel(cam, x, y, x_offset, y_offset);
            c = color_at(w, r, 5);
            sum[0] += c->arr[0];
            sum[1] += c->arr[1];
            sum[2] += c->arr[2];
            ray_free(r);
            color_free(c);
        }
    }

    c = color(sum[0] / total_steps, sum[1] / total_steps, sum[2] / total_steps);
    return c;
}

Canvas
render(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter)
{
    int i,j,k;

    Canvas image = canvas_alloc(cam->hsize, cam->vsize);

    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            Color c = pixel_multi_sample(cam, w, (double)i, (double)j, usteps, vsteps, jitter);
            //Color c = pixel_single_sample(cam, w, i, j);
            canvas_write_pixel(image, i, j, c);
            color_free(c);
        }
        k += 1;
        //printf("Wrote %d rows out of %lu\n", k, cam->vsize);
    }
    return image;
}

Color
color_at(World w, Ray r, size_t remaining)
{
    Intersections xs = intersect_world(w, r);
    Intersection i = hit(xs, false);
    if (i == NULL) {
        intersections_free(xs);
        return color(0,0,0);
    }

    Computations comps = prepare_computations(i, r, xs);
    Color c = shade_hit(w, comps, remaining);

    intersections_free(xs);

    return c;
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
    Shape itr;
    Intersections xs = intersections_empty(64);
    int i;

    for (itr = w->shapes, i = 0;
         i < w->shapes_num;
         itr++, i++) {
        Intersections xs_1 = itr->intersect(itr, r);
        if (xs_1->num == 0) {
            intersections_free(xs_1);
            continue;
        }

        // realloc and copy xs
        if (xs_1->num + xs->num >= xs->array_len) {
            Intersections xs_2 = intersections_empty(2 * xs->array_len);
            memcpy(xs_2->xs, xs->xs, xs->array_len * sizeof(struct intersection));
            xs_2->num = xs->num;
            intersections_free(xs);
            xs = xs_2;
        }

        // copy from xs_1 into xs + xs->num
        memcpy(xs->xs + xs->num, xs_1->xs, xs_1->num * sizeof(struct intersection));
        xs->num += xs_1->num;

        intersections_free(xs_1);
    }

    if (xs->num > 1) {
        intersections_sort(xs);
    }

    return xs;
}

Point
position_alloc(Ray ray, double t)
{
    Point position = point(ray->origin[0], ray->origin[1], ray->origin[2]);
    position->arr[0] += ray->direction[0] * t;
    position->arr[1] += ray->direction[1] * t;
    position->arr[2] += ray->direction[2] * t;

    return position;
}

Computations
computations(double t,
             Shape obj,
             Point p,
             Point over_point,
             Point under_point,
             Vector eyev,
             Vector normalv,
             Vector reflectv,
             bool inside)
{
    Computations comps = (Computations) malloc(sizeof(struct computations));
    comps->t = t;
    comps->obj = obj;
    comps->p = p;
    comps->eyev = eyev;
    comps->normalv = normalv;
    comps->reflectv = reflectv;
    comps->inside = inside;
    comps->over_point = over_point;
    comps->under_point = under_point;
    comps->n1 = 1.0;
    comps->n2 = 1.0;

    return comps;
}
 
Computations
prepare_computations(Intersection i, Ray r, Intersections xs)
{
    Point p = position_alloc(r, i->t);
    Vector n = i->object->normal_at(i->object, p, i);
    Vector eyev = vector(-r->direction[0], -r->direction[1], -r->direction[2]);
    Vector direction = vector(r->direction[0], r->direction[1], r->direction[2]);

    Point over_point = point_default();
    Point under_point = point_default();
    bool inside = false;

    Vector reflectv = vector_reflect_alloc(direction, n);

    if (vector_dot(n, eyev) < 0) {
        inside = true;
        vector_scale(n, -1);
    }

    over_point->arr[0] = p->arr[0] + n->arr[0] * EPSILON;
    over_point->arr[1] = p->arr[1] + n->arr[1] * EPSILON;
    over_point->arr[2] = p->arr[2] + n->arr[2] * EPSILON;

    under_point->arr[0] = p->arr[0] - n->arr[0] * EPSILON;
    under_point->arr[1] = p->arr[1] - n->arr[1] * EPSILON;
    under_point->arr[2] = p->arr[2] - n->arr[2] * EPSILON;

    Computations c = computations(i->t,
                                  i->object,
                                  p,
                                  over_point,
                                  under_point,
                                  eyev,
                                  n,
                                  reflectv,
                                  inside);

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
    vector_free(direction);

    return c;
}

void
computations_free(Computations comps)
{
    if (comps != NULL) {
        if (comps->p != NULL) {
            point_free(comps->p);
        }
        if (comps->over_point != NULL) {
            point_free(comps->over_point);
        }
        if (comps->under_point != NULL) {
            point_free(comps->under_point);
        }
        if (comps->eyev != NULL) {
            vector_free(comps->eyev);
        }
        if (comps->normalv != NULL) {
            vector_free(comps->normalv);
        }
        if (comps->reflectv != NULL) {
            vector_free(comps->reflectv);
        }
        free(comps);
    }
}

Color
reflected_color(World w, Computations comps, size_t remaining)
{
    if (remaining == 0 || equal(comps->obj->material->reflective, 0)) {
        return color(0,0,0);
    }
    Ray reflect_ray = ray_alloc(comps->over_point, comps->reflectv);
    Color c = color_at(w, reflect_ray, remaining - 1);
    color_scale(c, comps->obj->material->reflective);

    ray_free(reflect_ray);

    return c;
}

Color
refracted_color(World w, Computations comps, size_t remaining)
{
    if (remaining == 0 || equal(comps->obj->material->transparency,0)) {
        return color(0,0,0);
    }
    double n_ratio = comps->n1 / comps->n2;
    double cos_i = vector_dot(comps->eyev, comps->normalv);
    double sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);

    if (sin2_t > 1.0) {
        return color(0,0,0);
    }

    double cos_t = sqrt(1.0 - sin2_t);
    Vector t1 = vector(comps->normalv->arr[0], comps->normalv->arr[1], comps->normalv->arr[2]);
    vector_scale(t1, n_ratio * cos_i - cos_t);
    Vector t2 = vector(comps->eyev->arr[0], comps->eyev->arr[1], comps->eyev->arr[2]);
    vector_scale(t2, n_ratio);
    Vector direction = vector(t1->arr[0] - t2->arr[0],
                              t1->arr[1] - t2->arr[1],
                              t1->arr[2] - t2->arr[2]);

    Ray refracted_ray = ray_alloc(comps->under_point, direction);
    Color c = color_at(w, refracted_ray, remaining - 1);
    color_scale(c, comps->obj->material->transparency);

    ray_free(refracted_ray);
    vector_free(direction);
    vector_free(t2);
    vector_free(t1);

    return c;
}

double
schlick(Computations comps)
{
    double co = vector_dot(comps->eyev, comps->normalv);
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

Color
shade_hit(World w, Computations comps, size_t remaining)
{
    Color surface = color(0,0,0);
    Light itr;
    size_t i;

    for (i = 0, itr = w->lights; i < w->lights_num; i++, itr++) {
        double intensity = itr->intensity_at(itr, w, comps->over_point);
        Color c = lighting(comps->obj->material,
                           comps->obj,
                           itr,
                           comps->over_point,
                           comps->eyev,
                           comps->normalv,
                           intensity);

        color_accumulate(surface, c);
        color_free(c);
    }

    Color reflected = reflected_color(w, comps, remaining);
    Color refracted = refracted_color(w, comps, remaining);

    if (comps->obj->material->reflective > 0 && comps->obj->material->transparency > 0) {
        double reflectance = schlick(comps);
        color_scale(reflected, reflectance);
        color_scale(refracted, 1.0 - reflectance);
    }

    color_accumulate(surface, reflected);
    color_accumulate(surface, refracted);

    color_free(reflected);
    color_free(refracted);

    computations_free(comps);

    return surface;
}

Color
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity)
{
    Color pcolor;

    if (material->pattern != NULL) {
        pcolor = material->pattern->pattern_at_shape(material->pattern, shape, point);
    } else {
        pcolor = color_default();
        memcpy(pcolor->arr, material->color, sizeof(pcolor->arr));
    }
    Color effective_color = pcolor;

    pcolor->arr[0] *= light->intensity[0];
    pcolor->arr[1] *= light->intensity[1];
    pcolor->arr[2] *= light->intensity[2];

    Color ambient = color_default();
    memcpy(ambient->arr, effective_color->arr, sizeof(ambient->arr));
    color_scale(ambient, material->ambient);
    if (equal(shade_intensity,0.0)) {
        color_free(pcolor);
        return ambient;
    }

    Color acc = color(0,0,0);
    Color diffuse = color_default();
    int i;
    Points pts = light->light_surface_points(light);
    Point position;
    Vector diff = vector_default();
    Vector lightv = vector_default();
    Vector reflectv = vector_default();

    for (i = 0, position = pts->points; i < pts->points_num; i++, position++) {
        vector_from_points(position, point, diff);
        vector_normalize(diff, lightv);
        double light_dot_normal = vector_dot(lightv, normalv);
        if (light_dot_normal >= 0) {
            // diffuse
            diffuse->arr[0] = effective_color->arr[0];
            diffuse->arr[1] = effective_color->arr[1];
            diffuse->arr[2] = effective_color->arr[2];

            color_scale(diffuse, material->diffuse);
            color_scale(diffuse, light_dot_normal);
            color_accumulate(acc, diffuse);

            // specular
            vector_scale(lightv, -1);
            vector_reflect(lightv, normalv, reflectv);
            vector_scale(lightv, -1);

            double reflect_dot_eye = vector_dot(reflectv, eyev);
            if (reflect_dot_eye > 0) {
                double factor = pow(reflect_dot_eye, material->shininess);
                acc->arr[0] += light->intensity[0] * material->specular * factor;
                acc->arr[1] += light->intensity[1] * material->specular * factor;
                acc->arr[2] += light->intensity[2] * material->specular * factor;
            }
        }
    }

    double scaling = shade_intensity / (double)light->num_samples;
    color_scale(acc, scaling);
    color_accumulate(acc, ambient);

    vector_free(lightv);
    vector_free(diff);
    color_free(diffuse);
    color_free(ambient);
    color_free(pcolor);

    return acc;
}
