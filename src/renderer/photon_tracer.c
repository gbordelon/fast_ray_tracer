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
power_at(enum photon_map_type map_type, const World w, const Ray r, Color power, bool had_diffuse, bool had_specular, const size_t remaining);

// prepare_computations assumes the ray direction is from the camera but should be resusable

int
reflect_photon_diffuse(enum photon_map_type map_type, World w, Computations comps, bool had_specular, double average_diffuse_reflectance, size_t remaining)
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

    color_copy(diffuse_color, comps->over_Kd);
    diffuse_color[0] *= comps->photon_power[0] / average_diffuse_reflectance;
    diffuse_color[1] *= comps->photon_power[1] / average_diffuse_reflectance;
    diffuse_color[2] *= comps->photon_power[2] / average_diffuse_reflectance;
    // debug
    //color_copy(diffuse_color, comps->photon_power);

    ray_array(comps->over_point, random_diffuse, &reflect_ray);
    return power_at(map_type, w, &reflect_ray, diffuse_color, true, had_specular, remaining - 1);
}

int
reflect_photon_specular(enum photon_map_type map_type, World w, Computations comps, bool had_diffuse, double average_specular_reflectance, size_t remaining)
{
    if (!comps->obj->material->reflective) {
        return 0;
    }

    color_scale(comps->photon_power, 1.0 / average_specular_reflectance);

    struct ray reflect_ray;
    ray_array(comps->over_point, comps->reflectv, &reflect_ray);
    return power_at(map_type, w, &reflect_ray, comps->photon_power, had_diffuse, true, remaining - 1);
}

// TODO scale photon power by material->Tf triple but only once. Only as it either enters or leaves
// TODO how do i know? Can i use comps->inside ? 
int
refract_photon(enum photon_map_type map_type, World w, Computations comps, bool had_diffuse, double average_transmission_filter, size_t remaining)
{
    if (equal(comps->obj->material->Tr, 0)) {
        return 0;
    }

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

    color_scale(comps->photon_power, 1.0 / average_transmission_filter);

    return power_at(map_type, w, &refracted_ray, comps->photon_power, had_diffuse, true, remaining - 1);
}

int
photon_hit(enum photon_map_type map_type, World w, Computations comps, bool had_diffuse, bool had_specular, size_t remaining)
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

    Color diffuse_color, specular_color, transmission_filter;
    color_copy(diffuse_color, comps->over_Kd);
    color_copy(specular_color, comps->over_Ks);
    color_copy(transmission_filter, comps->obj->material->Tf);

    // only store photons if the struck material is diffuse
    if (diffuse_color[0] > 0 || diffuse_color[1] > 0 || diffuse_color[2] > 0) {
        // only store photons with at least one specular reflection in the caustics map
        if (map_type == CAUSTIC) {
            if (had_specular) {
                pm_store(maps + 0, comps->photon_power, comps->p, comps->photon_ray.direction);
                return 1;
            }
        // never store the first diffuse hit in the global map
        } else if (map_type == GLOBAL) {
            if (had_diffuse) {
                pm_store(maps + 1, comps->photon_power, comps->p, comps->photon_ray.direction);
                hit += 1;
            }
        }
        // participating media map?
    }

    // russian roulette
    double r = drand48();
    double total_refract_reflect;
    double average_specular_reflectance = (specular_color[0] + specular_color[1] + specular_color[2]) / 3.0;
    double average_transmission_filter = (transmission_filter[0] + transmission_filter[1] + transmission_filter[2]) / 3.0;
    if (map_type == GLOBAL) {
        double average_diffuse_reflectance = (diffuse_color[0] + diffuse_color[1] + diffuse_color[2]) / 3.0;
        total_refract_reflect = average_diffuse_reflectance +
                                average_specular_reflectance +
                                average_transmission_filter;
        if (r * total_refract_reflect < average_diffuse_reflectance) {
            hit += reflect_photon_diffuse(map_type, w, comps, had_specular, average_diffuse_reflectance, remaining);
        } else if (r * total_refract_reflect < average_diffuse_reflectance + average_specular_reflectance) {
            hit += reflect_photon_specular(map_type, w, comps, had_diffuse, average_specular_reflectance, remaining);
        } else if (r * total_refract_reflect < average_diffuse_reflectance + average_specular_reflectance + average_transmission_filter) {
            hit += refract_photon(map_type, w, comps, had_diffuse, average_transmission_filter, remaining);
        }
    } else if (map_type == CAUSTIC) {
        total_refract_reflect = average_specular_reflectance +
                                average_transmission_filter;
        if (r * total_refract_reflect < average_specular_reflectance) {
            hit += reflect_photon_specular(map_type, w, comps, had_diffuse, average_specular_reflectance, remaining);
        } else if (r * total_refract_reflect < average_specular_reflectance + average_transmission_filter) {
            hit += refract_photon(map_type, w, comps, had_diffuse, average_transmission_filter, remaining);
        }
    }
    // else participating media

    return hit;
}

int
power_at(enum photon_map_type map_type, const World w, const Ray r, Color power, bool had_diffuse, bool had_specular, const size_t remaining)
{
    Intersections xs = intersect_world(w, r, true); // TODO this may be wrong
    Intersection i = hit(xs, true);
    struct computations comps;
    int hit = 0;

    if (i != NULL) {
        prepare_computations(i, r, power, xs, &comps, &(w->container));
        hit = photon_hit(map_type, w, &comps, had_diffuse, had_specular, remaining);
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
            hit = power_at(CAUSTIC, w, &r, itr->intensity, false, false, 5);
            j += hit;
        }

        global_total = 0;
        for (j = 0; j < itr->num_photons;) {
            itr->emit_photon(itr, &r); // get ray from light
            hit = power_at(GLOBAL, w, &r, itr->intensity, false, false, 5);
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
