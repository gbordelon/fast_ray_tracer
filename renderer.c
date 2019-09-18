#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "renderer.h"

#include "linalg.h"
#include "shapes.h"
#include "sphere.h"
#include "plane.h"
#include "cube.h"

Color
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity);


/*
 *  set the cache fields when alloc'ing the first time
 *  no need to alloc after the first call for point light
 */
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
point_light_intensity_at(Light light, World w, Point p)
{
    if (is_shadowed(w, light->u.point.position, p)) {
        return 0.0;
    }
    return 1.0;
}

Light
point_light(Point p, Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    // null check l
    l->type = POINT_LIGHT;
    memcpy(l->intensity, intensity, sizeof(l->intensity));
    l->num_samples = 1;
    memcpy(l->u.point.position, p, sizeof(l->u.point.position));
    l->light_surface_points = point_light_surface_points;
    l->intensity_at = point_light_intensity_at;


    // populate surface_points_cache
    l->surface_points_cache = NULL;
    l->surface_points_cache = point_light_surface_points(l);
    return l;
}

Camera
camera(size_t hsize,
       size_t vsize,
       double field_of_view,
       Matrix transform)
{
    Camera c = (Camera) malloc(sizeof(struct camera));
    // null check c

    c->hsize = hsize;
    c->vsize = vsize;
    c->field_of_view = field_of_view;


    double half_view = tan(field_of_view / 2.0);
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
    Color c = color(1, 1, 1);
    Light l = point_light(p, c);
    w->lights = l;
    w->lights_num = 1;
    Shape shapes = (Shape) malloc(3 * sizeof(struct shape));
    Shape s1 = shapes;
    Shape s2 = shapes + 1;
    Shape s3 = shapes + 2;
    cube(s1);
    sphere(s2);
    plane(s3);

    s1->material->color[0] = 0.8;
    s1->material->color[1] = 1.0;
    s1->material->color[2] = 0.6;
    s1->material->diffuse = 0.7;
    s1->material->specular = 0.2;
    //s1->material->reflective = 1.0;
    //s1->material->refractive_index = 1.0;
    //s1->material->transparency = 1.0;
    //s1->material->casts_shadow = false;

    Matrix trans1 = matrix_translate_alloc(0.0, 1.0, 2.0);
    Matrix scale = matrix_scale_alloc(0.5, 0.5, 0.5);
    Matrix rotate = matrix_rotate_x_alloc(-0.1);
    Matrix trans2 = matrix_translate_alloc(0.0, -3.0, 0.0);

    shape_set_transform(s2, matrix_multiply_alloc(trans1, scale));
    shape_set_transform(s3, matrix_multiply_alloc(trans2, rotate));

    w->shapes = shapes;
    w->shapes_num = 3;

    //matrix_free(scale);
    //matrix_free(trans1);

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
ray_for_pixel(Camera cam, size_t px, size_t py)
{
    double xoffset = ((double)px + 0.5) * cam->pixel_size;
    double yoffset = ((double)py + 0.5) * cam->pixel_size;
    double world_x = cam->half_width - xoffset;
    double world_y = cam->half_height - yoffset;
    Matrix inv = cam->transform_inverse;
    Point p = point(world_x, world_y, -1);

    Point pixel = matrix_point_multiply_alloc(inv, p);
    p->arr[0] = 0;
    p->arr[1] = 0;
    p->arr[2] = 0;
    Point origin = matrix_point_multiply_alloc(inv, p);
    Vector v = vector_from_points_alloc(pixel, origin);
    Vector direction = vector_normalize_alloc(v);

    Ray r = ray_alloc(origin, direction);

    point_free(p);
    vector_free(v);

    return r;
}

Canvas
render(Camera cam, World w)
{
    int i,j,k;

    Canvas image = canvas_alloc(cam->hsize, cam->vsize);

    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            Ray r = ray_for_pixel(cam, i, j);
            Color c = color_at(w, r, 5);
            canvas_write_pixel(image, i, j, c);
            color_free(c);
            ray_free(r);
        }
        //printf("\n");
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
        return color(0,0,0);
    }
    Computations comps = prepare_computations(i, r, xs);

    intersections_free(xs);

    return shade_hit(w, comps, remaining);
}

int
sort_intersections(const void *p, const void *q)
{
    double l = ((Intersection)p)->t;
    double r = ((Intersection)q)->t;
    return round(l - r);
}

Intersections
intersect_world(World w, Ray r)
{
    Shape itr;
    Intersections xs = intersections_empty(512);
    int i;

    for (itr = w->shapes, i = 0;
         i < w->shapes_num;
         itr++, i++) {
        Intersections xs_1 = itr->intersect(itr, r);
        if (xs_1->num + xs->num >= xs->array_len) {
            // realloc and copy xs
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

    // sort xs by xs->xs->t ascending
    qsort((void*)xs->xs, xs->num, sizeof(struct intersection), sort_intersections);
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
    Vector neg_r_direction = vector(-r->direction[0], -r->direction[1], -r->direction[2]);
    Vector reflectv = vector(r->direction[0], r->direction[1], r->direction[2]);

    Point over_point = point_default();
    Point under_point = point_default();
    bool inside = false;

    Vector reflectv2 = vector_reflect_alloc(reflectv, n);

    if (vector_dot(n, neg_r_direction) < 0) {
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
                                  neg_r_direction,
                                  n,
                                  reflectv2,
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
    vector_free(reflectv);

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
    if (remaining == 0 || comps->obj->material->reflective == 0) {
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
    if (remaining == 0 || comps->obj->material->transparency == 0) {
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
        color_free(c);
    }
    computations_free(comps);

    return surface;
}

Color
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity)
{
    Color pcolor = color_default();
    Color effective_color = pcolor;

    if (material->pattern != NULL) {
        // nothing for now
        // pcolor = material.pattern.pattern_at_shape(shape, point)
    } else {
        memcpy(pcolor->arr, material->color, sizeof(pcolor->arr));
    }

    pcolor->arr[0] *= light->intensity[0];
    pcolor->arr[1] *= light->intensity[1];
    pcolor->arr[2] *= light->intensity[2];

    Color ambient = color_default();
    memcpy(ambient->arr, effective_color->arr, sizeof(ambient->arr));
    color_scale(ambient, material->ambient);
    if (shade_intensity == 0.0) {
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
