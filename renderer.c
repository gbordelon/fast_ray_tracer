#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "linalg.h"
#include "renderer.h"

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

    c->pixel_size = c->half_width * 2.0 / (double)hsize;
    c->transform = transform;

    return c;
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
    vector_free(v);
    vector_free(forward);
    vector_free(upn);
    vector_free(left);
    vector_free(true_up);
    matrix_free(orientation);
    matrix_free(m);
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
    Point p = point(-1, 10, -10);
    Color c = color(1, 1, 1);
    Light l = point_light(p, c);
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
    s1->material->transparency = 1.0;

    matrix_scale(0.5, 0.5, 0.5, s2->transform);

    w->shapes = shapes;
    w->shapes_num = 2;

    return w;
}

bool
is_shadowed(World w, double light_position[4], Point pt)
{
    Vector v = vector_from_arrays_alloc(light_position, pt->arr);
    double distance = vector_magnitude(v);
    Vector direction = vector_normalize_alloc(v);
    Ray r = ray_alloc(pt, direction);
    Intersections xs = intersect_world(w, r, true);
    Intersection h = hit(xs);
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
    // TODO cache the inverse
    Matrix inv = matrix_inverse_alloc(cam->transform);
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
    matrix_free(inv);

    return r;
}

Canvas
render(Camera cam, World w)
{
    int i,j,k;

    Canvas image = canvas_alloc(cam->hsize, cam->vsize);
/*
    >>> w = default_world()
    >>> r = ray(point(0,0,-5), vector(0,1,0))
    >>> c = color_at(w,r)
    >>> np.isclose(c, color(0,0,0))
    array([ True,  True,  True])

    >>> w = default_world()
    >>> r = ray(point(0,0,-5), vector(0,0,1))
    >>> c = color_at(w,r)
    >>> np.isclose(c, color(0.38066, 0.47583, 0.28549589))
    array([ True,  True,  True])

    >>> w = default_world()
    >>> outer = w.contains[0]
    >>> outer.material.ambient = 1.0
    >>> inner = w.contains[1]
    >>> inner.material.ambient = 1.0
    >>> r = ray(point(0,0,0.75), vector(0,0,-1))
    >>> c = color_at(w,r)
    >>> c == inner.material.color
    array([ True,  True,  True])
*/
    Point pt = point(0.0, 0.0, -5.0);
    Vector v = vector( 0.0, 0.0, 1.0);
    Ray r2 = ray_alloc(pt, v);
    Color c = color_at(w, r2, 5);
    char buf[256];
    color_to_string(buf, 256, c);


    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            Ray r = ray_for_pixel(cam, i, j);
            Color c = color_at(w, r, 5);
            canvas_write_pixel(image, i, j, c);
            color_free(c);
            ray_free(r);
        }
        k += 1;
        //printf("Wrote %d rows out of %lu\n", k, cam->vsize);
    }
    return image;
}

Color
color_at(World w, Ray r, size_t remaining)
{
    Intersections xs = intersect_world(w, r, false);
    Intersection i = hit(xs);
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
intersect_world(World w, Ray r, bool filter_shadow_casters)
{
    Shape itr;
    Intersections xs = intersections_empty(1024);
    size_t max_xs = 1024;
    int i;

    for (itr = w->shapes, i = 0;
         i < w->shapes_num;
         itr++, i++) {
        if (!filter_shadow_casters && itr->material->casts_shadow) {
            Intersections xs_1 = itr->intersect(itr, r);
            if (xs_1->num + xs->num >= max_xs) {
                // realloc and copy xs
                max_xs *= 2;
                Intersections xs_2 = intersections_empty(max_xs);
                memcpy(xs_2->xs, xs->xs, xs->array_len * sizeof(struct intersection));
                xs_2->array_len = max_xs;
                xs_2->num = xs->num;
                intersections_free(xs);
                xs = xs_2;
            }
            // copy from xs_1 into xs + xs->len
            memcpy(xs->xs + xs->num, xs_1->xs, xs_1->num * sizeof(struct intersection));
            xs->num += xs_1->num;

            intersections_free(xs_1);
        }
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
    comps->n1 = 0.0;
    comps->n2 = 0.0;

    return comps;
}
 
/*
    c.over_point = c.point + c.normalv * EPSILON
    c.under_point = c.point - c.normalv * EPSILON

    containers = []
    for i in xs:
        if i == intersection:
            if len(containers) == 0:
                c.n1 = 1.0
            else:
                c.n1 = containers[-1].material.refractive_index

        if i.object in containers:
            containers.remove(i.object)
        else:
            containers.append(i.object)

        if i == intersection:
            if len(containers) == 0:
                c.n2 = 1.0
            else:
                c.n2 = containers[-1].material.refractive_index
            break

    return c
*/
Computations
prepare_computations(Intersection i, Ray r, Intersections xs)
{
    Point p = position_alloc(r, i->t);
    Vector n = i->object->normal_at(i->object, p, i);
    Vector neg_r_direction = vector(-r->direction[0], -r->direction[1], -r->direction[2]);
    Vector reflectv = vector(r->direction[0], r->direction[1], r->direction[2]);

    Point over_point = point_default();
    Point under_point = point_default();

    over_point->arr[0] = p->arr[0] + n->arr[0] * EPSILON;
    over_point->arr[1] = p->arr[0] + n->arr[1] * EPSILON;
    over_point->arr[2] = p->arr[0] + n->arr[2] * EPSILON;

    under_point->arr[0] = p->arr[0] - n->arr[0] * EPSILON;
    under_point->arr[1] = p->arr[0] - n->arr[1] * EPSILON;
    under_point->arr[2] = p->arr[0] - n->arr[2] * EPSILON;

    Computations c = computations(i->t,
                                  i->object,
                                  p,
                                  over_point,
                                  under_point,
                                  neg_r_direction,
                                  n,
                                  vector_reflect_alloc(reflectv, n),
                                  false);

    if (vector_dot(c->normalv, c->eyev) < 0) {
        c->inside = true;
        c->normalv->arr[0] = -c->normalv->arr[0];
        c->normalv->arr[1] = -c->normalv->arr[1];
        c->normalv->arr[2] = -c->normalv->arr[2];
    }

    // TODO containers for c->n1 and c->n2

    vector_free(reflectv);

    return c;
}

/*
    if remaining == 0 or comps.object.material.reflective == 0:
        return color(0,0,0)

    reflect_ray = shapes.ray(comps.over_point, comps.reflectv)
    c = color_at(world, reflect_ray, remaining - 1)

    return c * comps.object.material.reflective
*/
Color
reflected_color(World w, Computations comps, size_t remaining)
{
    return color(0,0,0);
}

/*

    if comps.object.material.transparency == 0 or remaining == 0:
        return color(0,0,0)

    n_ratio = comps.n1 / comps.n2
    cos_i = dot(comps.eyev, comps.normalv)
    sin2_t = n_ratio ** 2 * (1 - cos_i ** 2)

    if sin2_t > 1.0:
        return color(0,0,0)

    cos_t = np.sqrt(1.0 - sin2_t)
    direction = comps.normalv * (n_ratio * cos_i - cos_t) - \
                comps.eyev * n_ratio

    refracted_ray = shapes.ray(comps.under_point, direction)
    c = color_at(world, refracted_ray, remaining - 1) * \
        comps.object.material.transparency
    return c
*/
Color
refracted_color(World w, Computations comps, size_t remaining)
{
    return color(0,0,0);
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

        if (comps->obj->material->reflective > 0 && comps->obj->material) {
            double reflectance = schlick(comps);
            color_scale(reflected, reflectance);
            color_scale(refracted, 1.0 - reflectance);
            color_accumulate(surface, reflected);
            color_accumulate(surface, refracted);
        } else {
            color_accumulate(surface, reflected);
            color_accumulate(surface, refracted);
        }
        color_free(c);
        color_free(reflected);
        color_free(refracted);
    }

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
            diffuse->arr[0] = effective_color->arr[0];
            diffuse->arr[1] = effective_color->arr[1];
            diffuse->arr[2] = effective_color->arr[2];

            color_scale(diffuse, material->diffuse);
            color_scale(diffuse, light_dot_normal);
            color_accumulate(acc, diffuse);

            vector_scale(lightv, -1);
            vector_reflect(lightv, normalv, reflectv);
            vector_scale(lightv, -1);
            double reflect_dot_eye = vector_dot(reflectv, eyev);
            // specular
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
