#include <stdbool.h>
#include <stdlib.h>

#include "../libs/sampler/sampler.h"
#include "../libs/linalg/linalg.h"

#include "../color/color.h"
#include "../renderer/world.h"
#include "../renderer/renderer.h"

#include "light.h"

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

#define AREA_LIGHT_CACHE_SIZE 65536
Points
area_light_surface_points(Light light)
{
    if (light->surface_points_cache == NULL) {
        int u, v, i;
        double jitter[2];
        size_t index[2];
        Point point = {0.0, 0.0, 0.0, 1.0};
        Points pts = (Points) malloc(AREA_LIGHT_CACHE_SIZE * sizeof(struct pts));
        Points itr = pts;
        struct sampler sampler;

        sampler_2d(light->u.area.jitter, light->u.area.usteps, light->u.area.vsteps, sampler_default_constraint, &sampler);
        
        for (i = 0, itr = pts; i < AREA_LIGHT_CACHE_SIZE; i++, itr++) {
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
           Color intensity,
           double photon_count,
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
    l->photon_count = photon_count;


    l->light_surface_points = area_light_surface_points;
    l->intensity_at = area_light_intensity_at;

    // populate surface_points_cache
    l->surface_points_cache = NULL;
    area_light_surface_points(l);
}

Light
area_light_alloc(Point corner,
                 Vector full_uvec,
                 size_t usteps,
                 Vector full_vvec,
                 size_t vsteps,
                 bool jitter,
                 double photon_count,
                 Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    area_light(corner, full_uvec, usteps, full_vvec, vsteps, jitter, intensity, photon_count, l);
    return l;
}


void
point_light(Point p, Color intensity, double photon_count, Light l)
{
    l->type = POINT_LIGHT;
    color_copy(l->intensity, intensity);
    l->num_samples = 1;
    point_copy(l->u.point.position, p);
    l->light_surface_points = point_light_surface_points;
    l->intensity_at = point_light_intensity_at;
    l->photon_count = photon_count;


    // populate surface_points_cache
    l->surface_points_cache = NULL;
    point_light_surface_points(l);
}

Light
point_light_alloc(Point p, double photon_count, Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    // null check l
    point_light(p, intensity, photon_count, l);
    return l;
}

Light
array_of_lights(size_t num)
{
    return (Light)malloc(num * sizeof(struct light));
}
