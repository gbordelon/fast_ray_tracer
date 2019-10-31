#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../libs/linalg/linalg.h"
#include "../libs/photon_map/pm.h"
#include "../libs/sampler/sampler.h"
#include "../libs/thpool/thpool.h"
#include "../libs/core_select/core_select.h"

#include "../pattern/pattern.h"
#include "../material/material.h"
#include "../shapes/shapes.h"
#include "../shapes/sphere.h"
#include "../shapes/plane.h"
#include "../shapes/cube.h"
#include "../shapes/cone.h"
#include "../shapes/cylinder.h"
#include "../shapes/triangle.h"
#include "../shapes/csg.h"
#include "../shapes/group.h"
#include "../shapes/toroid.h"

#include "renderer.h"
#include "world.h"
#include "camera.h"
#include "ray.h"

void color_at(const World w, const Ray r, const size_t remaining, Color res);
void lighting(Computations comps, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, bool use_ambient, bool use_diffuse, bool use_specular_highlights, Color res);
void lighting_microfacet(Computations comps, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, bool use_ambient, bool use_diffuse, bool use_specular_highlights, Color res);
void lighting_gi(Computations comps, Shape shape, Point point, Vector eyev, Vector normalv, PhotonMap *maps, Color res);
void lighting_caustics(Computations comps, Shape shape, Point point, Vector eyev, Vector normalv, PhotonMap *maps, Color res);
void shade_hit_gi(World w, Computations comps, Color res);

// for easier reading
static bool visualize_photon_map;
static bool include_direct;
static bool include_ambient;
static bool include_diffuse;
static bool include_spec_highlight;
static bool include_specular;
static size_t path_length;
static bool use_gi;
static bool use_caustics;
static bool use_final_gather;
static size_t final_gather_usteps;
static size_t final_gather_vsteps;
static size_t irradiance_estimate_num;
static double irradiance_estimate_radius;
static double irradiance_estimate_cone_filter_k;

void
setup_config(World w)
{
    visualize_photon_map = w->global_config->illumination.debug_visualize_photon_map;
    include_direct = w->global_config->illumination.include_direct;
    include_ambient = w->global_config->illumination.di.include_ambient;
    include_diffuse = w->global_config->illumination.di.include_diffuse;
    include_spec_highlight = w->global_config->illumination.di.include_specular_highlight;
    include_specular = w->global_config->illumination.di.include_specular;
    use_gi = w->global_config->illumination.include_global || visualize_photon_map;
    use_caustics = w->global_config->illumination.gi.include_caustics;
    use_final_gather = w->global_config->illumination.gi.include_final_gather;
    final_gather_usteps = w->global_config->illumination.gi.usteps;
    final_gather_vsteps = w->global_config->illumination.gi.vsteps;
    irradiance_estimate_num = w->global_config->illumination.gi.irradiance_estimate_num;
    irradiance_estimate_radius = w->global_config->illumination.gi.irradiance_estimate_radius;
    irradiance_estimate_cone_filter_k =  w->global_config->illumination.gi.irradiance_estimate_cone_filter_k;
    path_length = w->global_config->illumination.di.path_length;
}

bool
is_shadowed(World w, Point light_position, Point pt)
{
    Vector v;
    vector_from_points(light_position, pt, v); // from pt to light

    double distance = vector_magnitude(v); // distance between pt and light

    Vector direction;
    vector_normalize(v, direction);

    struct ray r;
    ray_array(pt, direction, &r);

    Intersections xs = intersect_world(w, &r, true);

    Intersection h = hit(xs, true);
    bool retval = h != NULL && h->t < distance;

    return retval;
}

void
ray_for_pixel(Camera cam, size_t px, size_t py, double xy_jitter[2], Ray res)
{
    double xoffset = ((double)px + xy_jitter[0]) * cam->pixel_size;
    double yoffset = ((double)py + xy_jitter[1]) * cam->pixel_size;
    double world_x = cam->half_width - xoffset;
    double world_y = cam->half_height - yoffset;
    Matrix inv;
    Point origin;
    Vector direction;
    Point p;
    Point pixel;
    Vector v;

    matrix_copy(cam->transform_inverse, inv);

    p[0] = world_x;
    p[1] = world_y;
    p[2] = -cam->canvas_distance;
    p[3] = 1.0;
    matrix_point_multiply(inv, p, pixel);

    p[0] = 0;
    p[1] = 0;
    sample_aperture(p, px, py, &cam->aperture);

    p[0] *= cam->aperture.size;
    p[1] *= cam->aperture.size;
    p[2] = 0;
    matrix_point_multiply(inv, p, origin);
    vector_from_points(pixel, origin, v);
    vector_normalize(v, direction);

    ray_array(origin, direction, res);
}

void
pixel_multi_sample(Camera cam, World w, size_t x, size_t y, size_t usteps, size_t vsteps, Sampler sampler, ColorTriple res)
{
    double jitter[2];
    size_t u, v, index[2];
    double total_steps = (double)usteps * (double)vsteps;
    struct ray r;
    ColorTriple c;
    ColorTriple acc;

    color_default(ambient_from_triple(acc));
    color_default(diffuse_from_triple(acc));
    color_default(specular_from_triple(acc));

    sampler->reset(sampler);

    for (v = 0; v < vsteps; v++) {
        index[1] = v;
        for (u = 0; u < usteps; u++) {
            index[0] = u;
            sampler->get_point(sampler, index, jitter);
            r.origin[0] = 0;
            r.origin[1] = 0;
            r.origin[2] = 0;
            r.origin[3] = 1;
            r.direction[0] = 0;
            r.direction[1] = 0;
            r.direction[2] = 0;
            r.direction[0] = 0;

            ray_for_pixel(cam, x, y, jitter, &r);

            color_default(ambient_from_triple(c));
            color_default(diffuse_from_triple(c));
            color_default(specular_from_triple(c));
            color_at(w, &r, path_length, c);

            color_accumulate(ambient_from_triple(acc), ambient_from_triple(c));
            color_accumulate(diffuse_from_triple(acc), diffuse_from_triple(c));
            color_accumulate(specular_from_triple(acc), specular_from_triple(c));
        }
    }

    color_scale(ambient_from_triple(acc), 1.0 / total_steps);
    color_scale(diffuse_from_triple(acc), 1.0 / total_steps);
    color_scale(specular_from_triple(acc), 1.0 / total_steps);

    color_copy(ambient_from_triple(res), ambient_from_triple(acc));
    color_copy(diffuse_from_triple(res), diffuse_from_triple(acc));
    color_copy(specular_from_triple(res), specular_from_triple(acc));
}

struct render_args {
    size_t y_start;
    size_t y_end;
    size_t usteps;
    size_t vsteps;
    bool jitter;
    int core_id;
    Camera cam;
    Canvas image;
};

void
render_multi_helper(World w, void *args)
{
    Camera cam = ((struct render_args *)args)->cam;
    Canvas image = ((struct render_args *)args)->image;
    size_t y_start = ((struct render_args *)args)->y_start;
    size_t y_end = ((struct render_args *)args)->y_end;
    size_t usteps = ((struct render_args *)args)->usteps;
    size_t vsteps = ((struct render_args *)args)->vsteps;
    bool jitter = ((struct render_args *)args)->jitter;
    //int core_id = ((struct render_args *)args)->core_id;

    //stick_this_thread_to_core(core_id);

    ColorTriple c, pixel_color;

    struct sampler sampler;
    sampler_2d(jitter, usteps, vsteps, sampler_default_constraint, &sampler);

    Color *buf = (Color *)malloc(cam->hsize * sizeof(Color));

    int i, j, k;
    for (j = y_start, k=1; j < y_end; ++j, ++k) {
        for (i = 0; i < cam->hsize; ++i) {
            color_default(ambient_from_triple(c));
            color_default(diffuse_from_triple(c));
            color_default(specular_from_triple(c));
            color_default(ambient_from_triple(pixel_color));
            color_default(diffuse_from_triple(pixel_color));
            color_default(specular_from_triple(pixel_color));
//            if (200 <= i && i <= 200 && 232 <= j && j <= 232) { // debug
                pixel_multi_sample(cam, w, i, j, usteps, vsteps, &sampler, c);
                // aggregate colors
                color_accumulate(pixel_color, ambient_from_triple(c));
                color_accumulate(pixel_color, diffuse_from_triple(c));
                color_accumulate(pixel_color, specular_from_triple(c));
                color_scale(pixel_color, 1.0 / 3.0);
//            }
            // record the color
            color_copy(*(buf+i), pixel_color);
        }
        printf("Wrote row %lu\n", y_start);
        canvas_write_pixels(image, 0, j, buf, cam->hsize);
    }

    sampler_free(&sampler);
    free(buf);
}

Canvas
render_multi(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter)
{
    setup_config(w);

    int i;
    size_t num_threads = w->global_config->threading.num_threads;
    Canvas image = canvas_alloc(cam->hsize, cam->vsize, false, NULL);
    World worlds = (World) malloc(num_threads * sizeof(struct world));
    struct render_args *args_array = (struct render_args *)malloc(cam->vsize * sizeof(struct render_args));
    struct render_args *debug = args_array;
    for (i = 0; i < num_threads; ++i) {
        world_copy(w, worlds + i);
    }

    threadpool thpool = thpool_init(num_threads, worlds);

    for (i = 0; i < cam->vsize; ++i) {
        debug->cam = cam;
        debug->image = image;
        debug->y_start = i;
        debug->y_end = i + 1;
        debug->usteps = usteps;
        debug->vsteps = vsteps;
        debug->jitter = jitter;
        debug->core_id = i % 4; // 4 cores on this Mac
        thpool_add_work(thpool, render_multi_helper, debug);
        debug++;
    }

    thpool_wait(thpool);
    thpool_destroy(thpool);

    free(args_array);
    // TODO recursive world free for each world
    free(worlds);

    return image;
}

Canvas
render(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter)
{
    setup_config(w);

    int i,j,k;
    ColorTriple c;
    Color pixel_color;

    Canvas image = canvas_alloc(cam->hsize, cam->vsize, false, NULL);
    struct sampler sampler;
    sampler_2d(jitter, usteps, vsteps, sampler_default_constraint, &sampler);

    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            color_default(ambient_from_triple(c));
            color_default(diffuse_from_triple(c));
            color_default(specular_from_triple(c));
            color_default(pixel_color);
            pixel_multi_sample(cam, w, i, j, usteps, vsteps, &sampler, c);
            // aggregate colors
            color_accumulate(pixel_color, ambient_from_triple(c));
            color_accumulate(pixel_color, diffuse_from_triple(c));
            color_accumulate(pixel_color, specular_from_triple(c));
            color_scale(pixel_color, 1.0 / 3.0);
            // record the color
            canvas_write_pixel(image, i, j, pixel_color);
        }
        k += 1;
        printf("Wrote %d rows out of %lu\n", k, cam->vsize);
    }

    return image;
}

void
color_at_gi(const World w, const Ray r, Color res)
{
    Intersections xs = intersect_world(w, r, false);
    Intersection i = hit(xs, false);
    struct computations comps;
    Color c;
    color_default(c);

    if (i != NULL) {
        Color diffuse;
        Point p;
        ray_position(r, i->t, p);
        if (i->object->material->map_Kd != NULL) {
            i->object->material->map_Kd->pattern_at_shape(i->object->material->map_Kd, i->object, p, diffuse);
        } else {
            color_copy(diffuse, i->object->material->Kd);
        }

        if (diffuse[0] > 0 || diffuse[1] > 0 || diffuse[2] > 0) {
            prepare_computations(i, r, c, xs, &comps, &(w->container)); // can probably do something simpler here
            shade_hit_gi(w, &comps, c);
        }
    }

    color_copy(res, c);
}

void
color_at(const World w, const Ray r, const size_t remaining, ColorTriple res)
{
    Intersections xs = intersect_world(w, r, false);
    Intersection i = hit(xs, false);
    struct computations comps;
    ColorTriple c;
    color_default(ambient_from_triple(c));
    color_default(diffuse_from_triple(c));
    color_default(specular_from_triple(c));

    if (i != NULL) {
        prepare_computations(i, r, diffuse_from_triple(c), xs, &comps, &(w->container));
        shade_hit(w, &comps, remaining, c);
    }

    color_copy(ambient_from_triple(res), ambient_from_triple(c));
    color_copy(diffuse_from_triple(res), diffuse_from_triple(c));
    color_copy(specular_from_triple(res), specular_from_triple(c));
}

void
prepare_computations(Intersection i, Ray r, Color photon_power, Intersections xs, Computations res, struct container *container)
{
    ray_array(r->origin, r->direction, &(res->photon_ray));
    color_copy(res->photon_power, photon_power);

    res->t = i->t;

    res->obj = i->object;

    ray_position(r, i->t, res->p);

    res->obj->normal_at(res->obj, res->p, i, res->normalv);

    vector_copy(res->eyev, r->direction);
    vector_scale(res->eyev, -1.0);

    res->inside = false;

    if (vector_dot(res->normalv, res->eyev) < 0) {
        res->inside = true;
        vector_scale(res->normalv, -1);
    }
    vector_reflect(r->direction, res->normalv, res->reflectv);

    res->over_point[0] = res->p[0] + res->normalv[0] * EPSILON;
    res->over_point[1] = res->p[1] + res->normalv[1] * EPSILON;
    res->over_point[2] = res->p[2] + res->normalv[2] * EPSILON;
    res->over_point[3] = 1.0;

    res->under_point[0] = res->p[0] - res->normalv[0] * EPSILON;
    res->under_point[1] = res->p[1] - res->normalv[1] * EPSILON;
    res->under_point[2] = res->p[2] - res->normalv[2] * EPSILON;
    res->under_point[3] = 1.0;

    res->n1 = 1.0;
    res->n2 = 1.0;

    int j, k;
    Intersection x;
    size_t container_len = 0;

    // check size of container array
    if (xs->num >= container->size) {
        size_t new_container_size = xs->num > (2 * container->size) ? xs->num : (2 * container->size);
        container->shapes = realloc(container->shapes, new_container_size * sizeof(Shape *));
        container->size = new_container_size;
    }

    for (j = 0, x = xs->xs; j < xs->num; x++, j++) {
        if (x == i) {// address compare should be okay.
            if (container_len > 0) {
                res->n1 = container->shapes[container_len-1]->material->Ni;
            }
        }

        for (k = 0; k < container_len; k++) {
            if (container->shapes[k] == x->object) {
                break;
            }
        }

        if (k < container_len) {
            // shift everything left one slot, overwriting container[index_of_object] first
            --container_len;
            for (; k < container_len; k++) {
                container->shapes[k] = container->shapes[k+1];
            }
        } else {
            container->shapes[container_len] = x->object;
            ++container_len;
        }

        if (x == i) {
            if (container_len > 0) {
                res->n2 = container->shapes[container_len-1]->material->Ni;
            }
            break;
        }
    }

    if (res->obj->material->map_Ka != NULL) {
        res->obj->material->map_Ka->pattern_at_shape(res->obj->material->map_Ka, res->obj, res->over_point, res->over_Ka);
    } else {
        color_copy(res->over_Ka, res->obj->material->Ka);
    }

    if (res->obj->material->map_Kd != NULL) {
        res->obj->material->map_Kd->pattern_at_shape(res->obj->material->map_Kd, res->obj, res->over_point, res->over_Kd);
    } else {
        color_copy(res->over_Kd, res->obj->material->Kd);
    }

    if (res->obj->material->map_Ks != NULL) {
        res->obj->material->map_Ks->pattern_at_shape(res->obj->material->map_Ks, res->obj, res->over_point, res->over_Ks);
    } else {
        color_copy(res->over_Ks, res->obj->material->Ks);
    }

    if (res->obj->material->map_refl != NULL) {
        res->obj->material->map_refl->pattern_at_shape(res->obj->material->map_refl, res->obj, res->over_point, res->over_refl);
    } else {
        color_copy(res->over_refl, res->obj->material->refl);
    }
/*
 * comps->reflectv is used for spec reflection not spec highlight
 * This conversion probably doesn't need to happen for spec highlight
𝑅 = 𝑟𝑒𝑓𝑙𝑒𝑐𝑡(𝑉, 𝑁);
𝑃 = 𝐼 + 𝑅 + 𝑠𝑐𝑎𝑙𝑒(𝑟𝑎𝑛𝑑𝑜𝑚𝐷𝑖𝑟𝑒𝑐𝑡𝑖𝑜𝑛𝐼𝑛𝐻𝑒𝑚𝑖𝑠𝑝ℎ𝑒𝑟𝑒𝐶𝑜𝑠𝑖𝑛𝑒𝑊𝑒𝑖𝑔ℎ𝑡𝑒𝑑(𝑅), 𝑠𝑝𝑒𝑐𝑢𝑙𝑎𝑟𝑖𝑡𝑦)
𝑅 = 𝑛𝑜𝑟𝑚𝑎𝑙𝑖𝑧𝑒(𝑃 − 𝐼);
*/

    if (res->obj->material->map_Ns != NULL) {
        Color tmp;
        res->obj->material->map_Ns->pattern_at_shape(res->obj->material->map_Ns, res->obj, res->over_point, tmp);
        res->over_Ns = tmp[0];
    } else {
        res->over_Ns = res->obj->material->Ns;
    }

    if (res->obj->material->map_d != NULL) {
        Color tmp;
        res->obj->material->map_d->pattern_at_shape(res->obj->material->map_d, res->obj, res->over_point, tmp);
        res->over_d = tmp[0]; // I hope the canvas reader stores this in the R channel...
    } else {
        res->over_d = 1.0 - res->obj->material->Tr;
    }
}

void
reflected_color(World w, Computations comps, size_t remaining, Color res)
{
    if (remaining == 0 || !comps->obj->material->reflective) {
        color_default(ambient_from_triple(res));
        color_default(diffuse_from_triple(res));
        color_default(specular_from_triple(res));
        return;
    }

    ColorTriple c;
    color_default(ambient_from_triple(c));
    color_default(diffuse_from_triple(c));
    color_default(specular_from_triple(c));

    struct ray reflect_ray;
    ray_array(comps->over_point, comps->reflectv, &reflect_ray);

    color_at(w, &reflect_ray, remaining - 1, c);

    ambient_from_triple(c)[0] *= comps->over_refl[0];
    ambient_from_triple(c)[1] *= comps->over_refl[1];
    ambient_from_triple(c)[2] *= comps->over_refl[2];

    diffuse_from_triple(c)[0] *= comps->over_refl[0];
    diffuse_from_triple(c)[1] *= comps->over_refl[1];
    diffuse_from_triple(c)[2] *= comps->over_refl[2];

    specular_from_triple(c)[0] *= comps->over_refl[0];
    specular_from_triple(c)[1] *= comps->over_refl[1];
    specular_from_triple(c)[2] *= comps->over_refl[2];

    color_accumulate(ambient_from_triple(res), ambient_from_triple(c));
    color_accumulate(diffuse_from_triple(res), diffuse_from_triple(c));
    color_accumulate(specular_from_triple(res), specular_from_triple(c));
}

void
refracted_color(World w, Computations comps, size_t remaining, ColorTriple res)
{
    if (remaining == 0 || comps->over_d <= 0.0) {
        color_default(ambient_from_triple(res));
        color_default(diffuse_from_triple(res));
        color_default(specular_from_triple(res));
        return;
    }

    double n_ratio = comps->n1 / comps->n2;
    double cos_i = vector_dot(comps->eyev, comps->normalv);
    double sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);

    if (sin2_t > 1.0) {
        color_default(ambient_from_triple(res));
        color_default(diffuse_from_triple(res));
        color_default(specular_from_triple(res));
        return;
    }

    ColorTriple c;
    color_default(ambient_from_triple(c));
    color_default(diffuse_from_triple(c));
    color_default(specular_from_triple(c));

    Vector t1;
    Vector t2;
    Vector direction;
    struct ray refracted_ray;

    double cos_t = sqrt(1.0 - sin2_t);
    vector_copy(t1, comps->normalv);
    vector_scale(t1, n_ratio * cos_i - cos_t);
    vector_copy(t2, comps->eyev);
    vector_scale(t2, n_ratio);

    vector(t1[0] - t2[0], t1[1] - t2[1], t1[2] - t2[2], direction);
    ray_array(comps->under_point, direction, &refracted_ray);

    color_at(w, &refracted_ray, remaining - 1, c);

    // Transmission filter
    ambient_from_triple(c)[0] *= comps->obj->material->Tf[0];
    ambient_from_triple(c)[1] *= comps->obj->material->Tf[1];
    ambient_from_triple(c)[2] *= comps->obj->material->Tf[2];

    diffuse_from_triple(c)[0] *= comps->obj->material->Tf[0];
    diffuse_from_triple(c)[1] *= comps->obj->material->Tf[1];
    diffuse_from_triple(c)[2] *= comps->obj->material->Tf[2];

    specular_from_triple(c)[0] *= comps->obj->material->Tf[0];
    specular_from_triple(c)[1] *= comps->obj->material->Tf[1];
    specular_from_triple(c)[2] *= comps->obj->material->Tf[2];

    // Dissolve
    ambient_from_triple(c)[0] *= comps->over_d;
    ambient_from_triple(c)[1] *= comps->over_d;
    ambient_from_triple(c)[2] *= comps->over_d;

    diffuse_from_triple(c)[0] *= comps->over_d;
    diffuse_from_triple(c)[1] *= comps->over_d;
    diffuse_from_triple(c)[2] *= comps->over_d;

    specular_from_triple(c)[0] *= comps->over_d;
    specular_from_triple(c)[1] *= comps->over_d;
    specular_from_triple(c)[2] *= comps->over_d;

    color_accumulate(ambient_from_triple(res), ambient_from_triple(c));
    color_accumulate(diffuse_from_triple(res), diffuse_from_triple(c));
    color_accumulate(specular_from_triple(res), specular_from_triple(c));
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

    double r0 = (comps->n1 - comps->n2) / (comps->n1 + comps->n2);
    r0 = r0 * r0;

    return r0 + (1.0 - r0) * (1.0 - co) * (1.0 - co) * (1.0 - co) * (1.0 - co) * (1.0 - co);
}

void
shade_hit_gi(World w, Computations comps, Color res)
{
    Color indirect;

    color_default(indirect);

    /* TODO allow for reflection and refraction */

    lighting_gi(comps,
                comps->obj,
                comps->over_point,
                comps->eyev,
                comps->normalv,
                w->photon_maps,
                indirect);

    color_copy(res, indirect);
    color_scale(res, M_PI);
}

void
final_gather(World w, Computations comps, Color res)
{
    Color total_power, c;
    struct ray r;
    double pdf_inv = 2 * M_PI; // initialize to a uniform hemisphere prob. dist. func.
    size_t num_rays = final_gather_usteps * final_gather_vsteps;
    size_t index[2];
    double rands[2];
    size_t num_samples = 0;
    struct sampler sampler;

    sampler_2d(true, final_gather_usteps, final_gather_vsteps, sampler_default_constraint, &sampler);
       
    point_copy(r.origin, comps->over_point);
    color_default(total_power);

    int i, j;
    for (j = 0; j < final_gather_vsteps; ++j) {
        index[1] = j;
        for (i = 0; i < final_gather_usteps; ++i) {
            index[0] = i;
            sampler.get_vector_hemisphere(&sampler, comps->normalv, true, index, rands, r.direction);
            color_default(c);
            color_at_gi(w, &r, c);

            color_scale(c, rands[0]); // scale by theta

            // accumulate the total power
            color_accumulate(total_power, c);
        }
    }

    color_scale(total_power, pdf_inv / (double)num_rays);

    res[0] = total_power[0] * comps->over_Kd[0];
    res[1] = total_power[1] * comps->over_Kd[1];
    res[2] = total_power[2] * comps->over_Kd[2];

    sampler_free(&sampler);
}

void
shade_hit(World w, Computations comps, size_t remaining, ColorTriple res)
{
    Light light;
    size_t i;
    ColorTriple c;
    ColorTriple surface;
    Color indirect;
    double intensity = 0;

    color_default(ambient_from_triple(surface));
    color_default(diffuse_from_triple(surface));
    color_default(specular_from_triple(surface));

    // shade the hit with direct ambient and direct specular highlights and maybe direct diffuse
    if (include_direct) {
        for (i = 0, light = w->lights; i < w->lights_num; i++, light++) {
            color_default(ambient_from_triple(c));
            color_default(diffuse_from_triple(c));
            color_default(specular_from_triple(c));
            intensity = light->intensity_at(light, w, comps->over_point);
            lighting_microfacet(comps,
                     comps->obj,
                     light,
                     comps->over_point,
                     comps->eyev,
                     comps->normalv,
                     intensity,
                     include_ambient,
                     include_diffuse,
                     include_spec_highlight,
                     c);
            color_accumulate(ambient_from_triple(surface), ambient_from_triple(c));
            color_accumulate(diffuse_from_triple(surface), diffuse_from_triple(c));
            color_accumulate(specular_from_triple(surface), specular_from_triple(c));
        }
    }

    // shade the hit with soft indirect diffuse and caustics
    if (use_gi && (comps->over_Kd[0] > 0 || comps->over_Kd[1] > 0 || comps->over_Kd[2] > 0)) {
        Color fgather;
        Color caustics;
        color_default(fgather);
        color_default(caustics);
        color_default(indirect);

        // approximate for the diffuse component, really only useful for visualizing the global pmap
        if (visualize_photon_map) {
            lighting_gi(comps,
                        comps->obj,
                        comps->over_point,
                        comps->eyev,
                        comps->normalv,
                        w->photon_maps,
                        indirect);
        }

        // final gather for soft indirect diffuse component
        if (use_final_gather) {
            final_gather(w, comps, fgather);
        }

        if (use_caustics) {
            lighting_caustics(comps,
                        comps->obj,
                        comps->over_point,
                        comps->eyev,
                        comps->normalv,
                        w->photon_maps,
                        caustics);
        }
        color_accumulate(ambient_from_triple(surface), indirect);
        color_accumulate(ambient_from_triple(surface), fgather);
        color_accumulate(ambient_from_triple(surface), caustics);

        // clamping
        double color_length = ambient_from_triple(surface)[0] + ambient_from_triple(surface)[1] + ambient_from_triple(surface)[2];
        if (color_length > 1.7320508075688772) { // sqrt(3) as max magnitude
            color_scale(ambient_from_triple(surface), 1.0 / color_length);
            color_scale(ambient_from_triple(surface), 1.7320508075688772);
        }
    }

    // shade the hit with specular reflections and refractions
    if (include_specular) {
        ColorTriple reflected;
        color_default(ambient_from_triple(reflected));
        color_default(diffuse_from_triple(reflected));
        color_default(specular_from_triple(reflected));

        reflected_color(w, comps, remaining, reflected);

        ColorTriple refracted;
        color_default(ambient_from_triple(refracted));
        color_default(diffuse_from_triple(refracted));
        color_default(specular_from_triple(refracted));

        refracted_color(w, comps, remaining, refracted);

        if (comps->obj->material->reflective && comps->over_d < 1.0) {
            double reflectance = schlick(comps);
            color_scale(ambient_from_triple(reflected), reflectance);
            color_scale(diffuse_from_triple(reflected), reflectance);
            color_scale(specular_from_triple(reflected), reflectance);

            color_scale(ambient_from_triple(refracted), 1.0 - reflectance);
            color_scale(diffuse_from_triple(refracted), 1.0 - reflectance);
            color_scale(specular_from_triple(refracted), 1.0 - reflectance);
        }

        color_accumulate(ambient_from_triple(surface), ambient_from_triple(reflected));
        color_accumulate(diffuse_from_triple(surface), diffuse_from_triple(reflected));
        color_accumulate(specular_from_triple(surface), specular_from_triple(reflected));


        if (comps->obj->material->Tr > 0.0 && comps->over_d > 0.0) {
            // Dissolve
            ambient_from_triple(surface)[0] *= 1.0 - comps->over_d;
            ambient_from_triple(surface)[1] *= 1.0 - comps->over_d;
            ambient_from_triple(surface)[2] *= 1.0 - comps->over_d;

            diffuse_from_triple(surface)[0] *= 1.0 - comps->over_d;
            diffuse_from_triple(surface)[1] *= 1.0 - comps->over_d;
            diffuse_from_triple(surface)[2] *= 1.0 - comps->over_d;

            specular_from_triple(surface)[0] *= 1.0 - comps->over_d;
            specular_from_triple(surface)[1] *= 1.0 - comps->over_d;
            specular_from_triple(surface)[2] *= 1.0 - comps->over_d;
        }

        color_accumulate(ambient_from_triple(surface), ambient_from_triple(refracted));
        color_accumulate(diffuse_from_triple(surface), diffuse_from_triple(refracted));
        color_accumulate(specular_from_triple(surface), specular_from_triple(refracted));
    }

    color_copy(ambient_from_triple(res), ambient_from_triple(surface));
    color_copy(diffuse_from_triple(res), diffuse_from_triple(surface));
    color_copy(specular_from_triple(res), specular_from_triple(surface));
}

void
lighting_caustics(Computations comps, Shape shape, Point point, Vector eyev, Vector normalv, PhotonMap *maps, Color res)
{
    Color intensity_estimate = {0.0, 0.0, 0.0, 0.0};
    Color caustic;

    color_copy(caustic, comps->over_Kd);


    if (!(caustic[0] > 0.0) && !(caustic[1] > 0.0) && !(caustic[2] > 0.0)) {
        return;
    }

    long num_photons_used = pm_irradiance_estimate(maps+0, intensity_estimate, point, eyev, irradiance_estimate_radius, irradiance_estimate_num, irradiance_estimate_cone_filter_k);

    if (num_photons_used > 0) {
        color_scale(intensity_estimate, 100.0 / ((double)num_photons_used));
    }

    if (visualize_photon_map) { // TODO investigate whether or not I should do this or always return the intensity_estimate
        color_copy(res, intensity_estimate);
    } else {
        double eye_dot_normal = vector_dot(eyev, normalv);

        caustic[0] *= intensity_estimate[0];
        caustic[1] *= intensity_estimate[1];
        caustic[2] *= intensity_estimate[2];

        color_scale(caustic, eye_dot_normal);
        color_accumulate(res, caustic);
    }
}

void
lighting_gi(Computations comps, Shape shape, Point point, Vector eyev, Vector normalv, PhotonMap *maps, Color res)
{
    Color intensity_estimate = {0.0, 0.0, 0.0, 0.0};
    Color diffuse;

    color_copy(diffuse, comps->over_Kd);

    if (!(diffuse[0] > 0.0) && !(diffuse[1] > 0.0) && !(diffuse[2] > 0.0)) {
        color_copy(res, diffuse);
        return;
    }

    long num_photons_used = pm_irradiance_estimate(maps+1, intensity_estimate, point, eyev, irradiance_estimate_radius, irradiance_estimate_num, irradiance_estimate_cone_filter_k);

    if (num_photons_used > 0) {
        color_scale(intensity_estimate, 10.0 * (double)irradiance_estimate_num / (double)num_photons_used);
    }

    double eye_dot_normal = vector_dot(eyev, normalv);
    if (visualize_photon_map) {
        color_copy(res, intensity_estimate);
    } else {
        diffuse[0] *= intensity_estimate[0];
        diffuse[1] *= intensity_estimate[1];
        diffuse[2] *= intensity_estimate[2];

        color_scale(diffuse, eye_dot_normal);
        color_copy(res, diffuse);
    }
}

void
lighting_microfacet(Computations comps, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, bool use_ambient, bool use_diffuse, bool use_specular_highlights, ColorTriple res)
{
    Color ambient;
    color_copy(ambient, comps->over_Ka);

    ambient[0] *= light->intensity[0];
    ambient[1] *= light->intensity[1];
    ambient[2] *= light->intensity[2];

    if (equal(shade_intensity, 0.0)) {
        if (use_ambient) {
            color_accumulate(ambient_from_triple(res), ambient);
        }
        return;
    }

    if (use_diffuse || use_specular_highlights) {
        Color c, diffuse, diffuse_acc, specular, specular_acc;
        int i;
        Vector diff, lightv, reflectv, halfwayv, tmp;
        Points pts = light->light_surface_points(light);
        double normal_dot_eye;

        if (use_diffuse) {
            color_copy(diffuse, comps->over_Kd);
        }
        if (use_specular_highlights) {
            color_copy(specular, comps->over_Ks);
            normal_dot_eye = vector_dot(normalv, eyev);
        }
        color_default(diffuse_acc);
        color_default(specular_acc);
        for (i = 0; i < pts->points_num; i++) {
            vector_from_points(*(pts->points + i), point, diff);
            vector_normalize(diff, lightv);
            double light_dot_normal = vector_dot(lightv, normalv);
            if (use_diffuse && light_dot_normal >= 0.0) {
                // diffuse
                color_copy(c, diffuse);
                c[0] *= light->intensity[0];
                c[1] *= light->intensity[1];
                c[2] *= light->intensity[2];

                color_scale(c, light_dot_normal);
                color_accumulate(diffuse_acc, c);
            }
            if (use_specular_highlights && light_dot_normal >= 0.0) {
                // specular
                double normal_dot_light = vector_dot(normalv, lightv);
                vector_default(tmp);
                tmp[0] = lightv[0] + eyev[0];
                tmp[1] = lightv[1] + eyev[1];
                tmp[2] = lightv[2] + eyev[2];
                tmp[3] = 0.0;
                vector_normalize(tmp, halfwayv);

                double normal_dot_halfway = fmax(0.0, vector_dot(normalv, halfwayv));
                double eye_dot_halfway_inv = 1.0 / fmax(0.0, vector_dot(eyev, halfwayv));
                double light_dot_halfway = vector_dot(lightv, halfwayv);

                double distribution_term = (comps->over_Ns + 2) * pow(normal_dot_halfway, comps->over_Ns) * 0.5 * M_1_PI;
                double geometry_coeff = 2.0 * normal_dot_halfway * eye_dot_halfway_inv;
                double geometry_term = fmin(1.0, fmin(geometry_coeff * normal_dot_eye,geometry_coeff * normal_dot_light));

                double factor = pow(1.0 - light_dot_halfway, 5.0);
                c[0] = specular[0] + (1.0 - specular[0]) * factor;
                c[1] = specular[1] + (1.0 - specular[1]) * factor;
                c[2] = specular[2] + (1.0 - specular[2]) * factor;
                double brdf_factor = distribution_term * geometry_term / (4.0 * normal_dot_light * normal_dot_eye);
                specular_acc[0] += c[0] * light->intensity[0] * brdf_factor;
                specular_acc[1] += c[1] * light->intensity[1] * brdf_factor;
                specular_acc[2] += c[2] * light->intensity[2] * brdf_factor;
            }
        }
        color_accumulate(diffuse_from_triple(res), diffuse_acc);
        color_accumulate(specular_from_triple(res), specular_acc);
        double scaling = shade_intensity / (double)light->num_samples;
        color_scale(diffuse_from_triple(res), scaling);
        color_scale(specular_from_triple(res), scaling);
    }

    if (use_ambient) {
        color_accumulate(ambient_from_triple(res), ambient);
    }
}

void
lighting(Computations comps, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, bool use_ambient, bool use_diffuse, bool use_specular_highlights, ColorTriple res)
{
    Color ambient;
    color_copy(ambient, comps->over_Ka);

    ambient[0] *= light->intensity[0];
    ambient[1] *= light->intensity[1];
    ambient[2] *= light->intensity[2];

    if (equal(shade_intensity, 0.0)) {
        if (use_ambient) {
            color_accumulate(ambient_from_triple(res), ambient);
        }
        return;
    }

    if (use_diffuse || use_specular_highlights) {
        bool phong_or_blinn_phong = false;
        Color diffuse, specular, diffuse_acc, specular_acc, c;
        int i;
        Vector diff, lightv, reflectv, halfwayv, tmp;
        Points pts = light->light_surface_points(light);

        color_default(diffuse_acc);
        color_default(specular_acc);

        if (use_diffuse) {
            color_copy(diffuse, comps->over_Kd);
        }
        if (use_specular_highlights) {
            color_copy(specular, comps->over_Ks);
        }
        for (i = 0; i < pts->points_num; i++) {
            vector_from_points(*(pts->points + i), point, diff);
            vector_normalize(diff, lightv);
            double light_dot_normal = vector_dot(lightv, normalv);
            if (light_dot_normal >= 0) {
                if (use_diffuse) {
                    // diffuse
                    color_copy(c, diffuse);
                    c[0] *= light->intensity[0];
                    c[1] *= light->intensity[1];
                    c[2] *= light->intensity[2];

                    color_scale(c, light_dot_normal);
                    color_accumulate(diffuse_acc, c);
                }
                if (use_specular_highlights) {
                    // specular
                    if (phong_or_blinn_phong) {
                        vector_scale(lightv, -1);
                        vector_reflect(lightv, normalv, reflectv);
                        vector_scale(lightv, -1);
                        double reflect_dot_eye = vector_dot(reflectv, eyev);
                        if (reflect_dot_eye > 0) {
                            double factor = pow(reflect_dot_eye, comps->obj->material->Ns);
                            specular_acc[0] += light->intensity[0] * specular[0] * factor;
                            specular_acc[1] += light->intensity[1] * specular[1] * factor;
                            specular_acc[2] += light->intensity[2] * specular[2] * factor;
                        }
                    } else {
                        tmp[0] = lightv[0] + eyev[0];
                        tmp[1] = lightv[1] + eyev[1];
                        tmp[2] = lightv[2] + eyev[2];
                        tmp[3] = 0.0;
                        vector_normalize(tmp, halfwayv);
                        double halfway_angle = vector_dot(normalv, halfwayv);
                        if (halfway_angle > 0) {
                            double factor = pow(halfway_angle, 4.0 * comps->over_Ns);
                            specular_acc[0] += light->intensity[0] * specular[0] * factor;
                            specular_acc[1] += light->intensity[1] * specular[1] * factor;
                            specular_acc[2] += light->intensity[2] * specular[2] * factor;
                        }
                    }
                }
            }
        }
        color_accumulate(diffuse_from_triple(res), diffuse_acc);
        color_accumulate(specular_from_triple(res), specular_acc);
        double scaling = shade_intensity / (double)light->num_samples;
        color_scale(diffuse_from_triple(res), scaling);
        color_scale(specular_from_triple(res), scaling);
    }

    if (use_ambient) {
        color_accumulate(ambient_from_triple(res), ambient);
    }
}
