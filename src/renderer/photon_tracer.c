#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "../libs/photon_map/pm.h"
#include "../libs/sampler/sampler.h"
#include "../color/color.h"
#include "../shapes/shapes.h"
#include "../pattern/pattern.h"

#include "renderer.h"
#include "world.h"
#include "photon_tracer.h"

enum reflect_type {
    NONE,
    DIFFUSE,
    SPECULAR
};

enum photon_map_type {
    CAUSTIC,
    GLOBAL,
    MEDIA
};

int
power_at(enum photon_map_type map_type, const World w, const Ray r, Color power, enum reflect_type type, const size_t remaining, struct container *container);

// prepare_computations assumes the ray direction is from the camera but should be resusable

int
reflect_photon_diffuse(enum photon_map_type map_type, World w, Computations comps, size_t remaining, struct container *container)
{
    Vector nt, nb, random_diffuse, sample;
    create_coordinate_system(comps->normalv, nt, nb);
    double r1 = drand48();
    double r2 = drand48();
    cosine_weighted_sample_hemisphere(r1, r2, sample);
    random_diffuse[0] = sample[0] * nb[0] + sample[1] * comps->normalv[0] + sample[2] * nt[0];
    random_diffuse[1] = sample[0] * nb[1] + sample[1] * comps->normalv[1] + sample[2] * nt[1];
    random_diffuse[2] = sample[0] * nb[2] + sample[1] * comps->normalv[2] + sample[2] * nt[2];
    random_diffuse[3] = 0;

    struct ray reflect_ray;
    Color diffuse_color;

    if (comps->obj->material->pattern != NULL) {
        comps->obj->material->pattern->pattern_at_shape(
                comps->obj->material->pattern,
                comps->obj,
                comps->p,
                diffuse_color);
    } else {
        color_copy(diffuse_color, comps->obj->material->color);
    }
    diffuse_color[0] *= comps->photon_power[0];
    diffuse_color[1] *= comps->photon_power[1];
    diffuse_color[2] *= comps->photon_power[2];

    ray_array(comps->over_point, random_diffuse, &reflect_ray);
    return power_at(map_type, w, &reflect_ray, diffuse_color, DIFFUSE, remaining - 1, container);
}

int
reflect_photon_specular(enum photon_map_type map_type, World w, Computations comps, size_t remaining, struct container *container)
{
    struct ray reflect_ray;
    ray_array(comps->over_point, comps->reflectv, &reflect_ray);
    return power_at(map_type, w, &reflect_ray, comps->photon_power, SPECULAR, remaining - 1, container);
}

int
refract_photon(enum photon_map_type map_type, World w, Computations comps, size_t remaining, struct container *container)
{
    double n_ratio = comps->n1 / comps->n2;
    double cos_i = vector_dot(comps->eyev, comps->normalv);
    double sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);

    if (sin2_t > 1.0) {
        return 0;
    }

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
    return power_at(map_type, w, &refracted_ray, comps->photon_power, SPECULAR, remaining - 1, container);
}

int
photon_hit(enum photon_map_type map_type, World w, Computations comps, enum reflect_type type, size_t remaining, struct container *container)
{
    PhotonMap *maps = w->photon_maps;
    int hit = 0;
    
    if (remaining <= 0) {
        return 0;
    }

    // don't bother recursing with shadow or dead photons
    if (comps->photon_power[0] <= 0 && comps->photon_power[1] <= 0 && comps->photon_power[2] <= 0) {
        return 0;
    }

    if (comps->obj->material->diffuse > 0) {
        if (map_type == CAUSTIC) {
            if (type == SPECULAR) { // only store photons with at least one specular reflection in the caustics map
                pm_store(maps + 0, comps->photon_power, comps->p, comps->photon_ray.direction);
                return 1;
            } else {
                return 0;
            }
        } else if (map_type == GLOBAL) {
            // never store the first diffuse hit in the global map
            if (type != NONE) {
                pm_store(maps + 1, comps->photon_power, comps->p, comps->photon_ray.direction);
                hit += 1;
            }
        }
    }

    // russian roulette
    double r = drand48();
    double total_refract_reflect = comps->obj->material->diffuse +
                                   comps->obj->material->reflective +
                                   comps->obj->material->transparency;

    if (map_type != CAUSTIC && (r * total_refract_reflect < comps->obj->material->diffuse)) {
        hit += reflect_photon_diffuse(map_type, w, comps, remaining, container);
    } else if (r * total_refract_reflect < comps->obj->material->diffuse + comps->obj->material->reflective) {
        hit += reflect_photon_specular(map_type, w, comps, remaining, container);
    } else if (r * total_refract_reflect < comps->obj->material->diffuse + comps->obj->material->reflective + comps->obj->material->transparency) {
        hit += refract_photon(map_type, w, comps, remaining, container);
    }

    return hit;
}

int
power_at(enum photon_map_type map_type, const World w, const Ray r, Color power, enum reflect_type type, const size_t remaining, struct container *container)
{
    Intersections xs = intersect_world(w, r);
    Intersection i = hit(xs, true);
    struct computations comps;
    int hit = 0;

    if (i != NULL) {
        prepare_computations(i, r, power, xs, &comps, container);
        hit = photon_hit(map_type, w, &comps, type, remaining, container);
    }

    return hit;
}

void
trace_photons(const World w, size_t num_maps)
{
    double total_lightness = 0;
    Color lab;
    struct ray r; 
    Light itr;
    int i, j, hit, global_total;
    struct container container;
    container.shapes = NULL;
    container.size = 0;
    PhotonMap *maps = w->photon_maps;
    size_t num_photons = maps->max_photons; // assume all three maps have the same photon count


    // associate total scene lightness with total number of photons
    for (i = 0, itr = w->lights; i < w->lights_num; ++i, ++itr) {
        rgb_to_lab(itr->intensity, lab);
        total_lightness += lab[0];
        itr->num_photons = 0;
    }

    // apportion the photons for each light
    for (i = 0, itr = w->lights; i < w->lights_num; ++i, ++itr) {
        rgb_to_lab(itr->intensity, lab);
        itr->num_photons = num_photons * lab[0] / total_lightness;
    }

    for (i = 0, itr = w->lights; i < w->lights_num; ++i, ++itr) {
        for (j = 0; j < itr->num_photons;) {
            itr->emit_photon(itr, &r); // get ray from light
            hit = power_at(CAUSTIC, w, &r, itr->intensity, NONE, 5, &container);
            j += hit;
        }

        global_total = 0;
        for (j = 0; j < itr->num_photons;) {
            itr->emit_photon(itr, &r); // get ray from light
            hit = power_at(GLOBAL, w, &r, itr->intensity, NONE, 5, &container);
            j += hit;
            global_total += hit;
        }
        for (j = 0; j < num_maps; ++j) {
            pm_scale_photon_power(maps + i, 1.0/(double)itr->num_photons);
        }
    }

    for (i = 0; i < num_maps; ++i) {
        pm_balance(maps + i);
    }
}

PhotonMap *
array_of_photon_maps(size_t num)
{
    return (PhotonMap *)malloc(num * sizeof(PhotonMap));
}
