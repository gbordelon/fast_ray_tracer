#include <stdbool.h>
#include <stdlib.h>

#include "../libs/sampler/sampler.h"
#include "../libs/linalg/linalg.h"

#include "../color/color.h"
#include "../renderer/world.h"
#include "../renderer/renderer.h"

#include "light.h"

void
area_light_emit_photon(Light l, Ray res)
{
    Vector normal, nt, nb, sample, direction;
    double r1 = drand48();
    double r2 = drand48();

    Vector uvec;
    vector_copy(uvec, l->u.area.uvec);
    vector_scale(uvec, l->u.area.usteps);

    Vector vvec;
    vector_copy(vvec, l->u.area.vvec);
    vector_scale(vvec, l->u.area.vsteps);

    vector_scale(uvec, r1);
    vector_scale(vvec, r2);

    res->origin[0] = l->u.area.corner[0] + uvec[0] + vvec[0];
    res->origin[1] = l->u.area.corner[1] + uvec[1] + vvec[1];
    res->origin[2] = l->u.area.corner[2] + uvec[2] + vvec[2];
    res->origin[3] = 1;

    vector_cross(l->u.area.uvec, l->u.area.vvec, normal); // may be backward
    create_coordinate_system(normal, nt, nb);

    r1 = drand48();
    r2 = drand48();
    cosine_weighted_sample_hemisphere(r1, r2, sample);
//    uniform_sample_hemisphere(r1, r2, sample);
//    direction[0] = sample[0] * nb[0] + sample[1] * normal[0] + sample[2] * nt[0];
//    direction[1] = sample[0] * nb[1] + sample[1] * normal[1] + sample[2] * nt[1];
//    direction[2] = sample[0] * nb[2] + sample[1] * normal[2] + sample[2] * nt[2];
//    direction[3] = 0;
//    vector_normalize(direction, res->direction);
    sample[1] *= -1.0;
    vector_normalize(sample, res->direction);
    //printf("direction %f %f %f\n", res->direction[0], res->direction[1], res->direction[2]);
}

void
hemisphere_light_emit_photon(Light l, Ray res)
{
    Vector nt, nb, sample, direction;
    double r1 = drand48();
    double r2 = drand48();
    point_copy(res->origin, l->u.hemi.position);

    create_coordinate_system(l->u.hemi.normal, nt, nb);

    r1 = drand48();
    r2 = drand48();
    cosine_weighted_sample_hemisphere(r1, r2, sample);
//    uniform_sample_hemisphere(r1, r2, sample);
    direction[0] = sample[0] * nb[0] + sample[1] * l->u.hemi.normal[0] + sample[2] * nt[0];
    direction[1] = sample[0] * nb[1] + sample[1] * l->u.hemi.normal[1] + sample[2] * nt[1];
    direction[2] = sample[0] * nb[2] + sample[1] * l->u.hemi.normal[2] + sample[2] * nt[2];
    direction[3] = 0;
    vector_normalize(direction, res->direction);
}

void
point_light_emit_photon(Light l, Ray res)
{
    Vector direction;
    direction[3] = 0;

    do {
        direction[0] = 2 * drand48() - 1;
        direction[1] = 2 * drand48() - 1;
        direction[2] = 2 * drand48() - 1;
    } while (direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2] > 1);

    point_copy(res->origin, l->u.point.position);
    vector_copy(res->direction, direction);
}

void
area_light_point_on_light(Light l, double uv_jitter[2], Point retval)
{
    Vector uvec;
    vector_copy(uvec, l->u.area.uvec);

    Vector vvec;
    vector_copy(vvec, l->u.area.vvec);

    vector_scale(uvec, uv_jitter[0]);
    vector_scale(vvec, uv_jitter[1]);

    retval[0] = l->u.area.corner[0] + uvec[0] + vvec[0];
    retval[1] = l->u.area.corner[1] + uvec[1] + vvec[1];
    retval[2] = l->u.area.corner[2] + uvec[2] + vvec[2];
}

void
construct_area_light_surface_points_cache(Light light, size_t cache_size)
{
    if (light->surface_points_cache == NULL) {
        int u, v, i;
        double jitter[2];
        size_t index[2];
        Point point = {0.0, 0.0, 0.0, 1.0};
        Points pts = (Points) malloc(cache_size * sizeof(struct pts));
        Points itr = pts;
        struct sampler sampler;

        sampler_2d(light->u.area.jitter, light->u.area.usteps, light->u.area.vsteps, sampler_default_constraint, &sampler);
        
        for (i = 0, itr = pts; i < cache_size; i++, itr++) {
            itr->points_num = light->num_samples;
            itr->points = (Point*) malloc(light->num_samples * sizeof(Point));

            sampler.reset(&sampler);

            for (v = 0; v < light->u.area.vsteps; ++v) {
                index[1] = v;
                for (u = 0; u < light->u.area.usteps; ++u) {
                    index[0] = u;
                    sampler.get_point(&sampler, index, jitter);
                    jitter[0] *= light->u.area.usteps;
                    jitter[1] *= light->u.area.vsteps;
                    
                    area_light_point_on_light(light, jitter, point);
                    point_copy(*(itr->points + v * light->u.area.usteps + u), point);
                }
            }
        }
        light->surface_points_cache = pts;
        light->surface_points_cache_len = cache_size;
    }
}

Points
area_light_surface_points(Light light)
{
    int choice = rand() % light->surface_points_cache_len;
    return light->surface_points_cache + choice;
}

Points
point_light_surface_points(Light light)
{
    if (light->surface_points_cache == NULL) {
        Points pts = (Points) malloc(sizeof(struct pts));
        pts->points_num = 1;
        pts->points = (Point*) malloc(pts->points_num * sizeof(Point));
        point_copy(*pts->points, light->u.point.position);
        light->surface_points_cache = pts;
        light->surface_points_cache_len = 1;
    }
    return light->surface_points_cache;
}

// TODO UV mapped sphere sampling like area light...
Points
hemisphere_light_surface_points(Light light)
{
    if (light->surface_points_cache == NULL) {
        Points pts = (Points) malloc(sizeof(struct pts));
        pts->points_num = 1;
        pts->points = (Point*) malloc(pts->points_num * sizeof(Point));
        point_copy(*pts->points, light->u.point.position);
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

    for (i = 0; i < pts->points_num; i++) {
        if (!is_shadowed(w, *(pts->points + i), p)) {
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
area_light(Point corner,
           Vector full_uvec,
           size_t usteps,
           Vector full_vvec,
           size_t vsteps,
           bool jitter,
           size_t cache_size,
           Color intensity,
           Light l)
{
    l->type = AREA_LIGHT;

    point_copy(l->u.area.corner, corner);

    vector_copy(l->u.area.uvec, full_uvec);
    vector_scale(l->u.area.uvec, 1.0 / (double) usteps);
    l->u.area.usteps = usteps;

    vector_copy(l->u.area.vvec, full_vvec);
    vector_scale(l->u.area.vvec, 1.0 / (double) vsteps);
    l->u.area.vsteps = vsteps;
    l->u.area.jitter = jitter;

    color_copy(l->intensity, intensity);
    l->num_samples = usteps * vsteps;


    l->light_surface_points = area_light_surface_points;
    l->intensity_at = area_light_intensity_at;
    l->emit_photon = area_light_emit_photon;

    // populate surface_points_cache
    l->surface_points_cache = NULL;
    construct_area_light_surface_points_cache(l, cache_size);
}

Light
area_light_alloc(Point corner,
                 Vector full_uvec,
                 size_t usteps,
                 Vector full_vvec,
                 size_t vsteps,
                 bool jitter,
                 size_t cache_size,
                 Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    area_light(corner, full_uvec, usteps, full_vvec, vsteps, jitter, cache_size, intensity, l);
    return l;
}


void
point_light(Point p, Color intensity, Light l)
{
    l->type = POINT_LIGHT;
    color_copy(l->intensity, intensity);
    l->num_samples = 1;
    point_copy(l->u.point.position, p);
    l->light_surface_points = point_light_surface_points;
    l->intensity_at = point_light_intensity_at;
    l->emit_photon = point_light_emit_photon;


    // populate surface_points_cache
    l->surface_points_cache = NULL;
    point_light_surface_points(l);
}

void
hemisphere_light(Point p, Vector normal, Color intensity, Light l)
{
    l->type = HEMISPHERE_LIGHT;
    color_copy(l->intensity, intensity);
    l->num_samples = 1;
    point_copy(l->u.hemi.position, p);
    vector_normalize(normal, l->u.hemi.normal);

    l->light_surface_points = hemisphere_light_surface_points;
    l->intensity_at = point_light_intensity_at; // TODO make this like area light
    l->emit_photon = hemisphere_light_emit_photon;


    // populate surface_points_cache
    l->surface_points_cache = NULL;
    hemisphere_light_surface_points(l);
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
