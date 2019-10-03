#include <stdbool.h>
#include <stdlib.h>

#include "../libs/sampler/sampler.h"
#include "../libs/linalg/linalg.h"

#include "../color/color.h"
#include "../renderer/world.h"
#include "../renderer/renderer.h"

#include "light.h"

void
area_light_point_on_light(Light l, size_t u, size_t v, double u_jitter, double v_jitter, Point retval)
{
    Vector uvec;
    vector_copy(uvec, l->u.area.uvec);

    Vector vvec;
    vector_copy(vvec, l->u.area.vvec);

    vector_scale(uvec, u_jitter);
    vector_scale(vvec, v_jitter);

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
        size_t n = light->u.area.vsteps;
        size_t m = light->u.area.usteps;
        size_t k;
        Point point = {0.0, 0.0, 0.0, 1.0};
        Points pts = (Points) malloc(AREA_LIGHT_CACHE_SIZE * sizeof(struct pts));
        Points itr = pts;
        static double *canonicalx = NULL;
        static double *canonicaly = NULL;

        if (canonicalx == NULL) {
            canonicalx = (double *) malloc(n * m * sizeof(double));
        }

        if (canonicaly == NULL) {
            canonicaly = (double *) malloc(n * m * sizeof(double));
        }

        for (i = 0, itr = pts; i < AREA_LIGHT_CACHE_SIZE; i++, itr++) {
            itr->points_num = light->num_samples;
            itr->points = (Point*) malloc(light->num_samples * sizeof(Point));

            // produce canonical representation
            for (v = 0; v < n; ++v) {
                for (u = 0; u < m; ++u) {
                    //canonicalx[(v * m + u)] = ((double)v + jitter_by(light->u.area.jitter)) / (double)n;
                    //canonicaly[(v * m + u)] = ((double)u + jitter_by(light->u.area.jitter)) / (double)m;
                    canonicalx[(v * m + u)] = ((double)u + ((double)v + jitter_by(light->u.area.jitter)) / (double)n) / (double)m;
                    canonicaly[(v * m + u)] = ((double)v + ((double)u + jitter_by(light->u.area.jitter)) / (double)m) / (double)n;

                    //canonicalx[(v * m + u)] = jitter_by(light->u.area.jitter);
                    //canonicaly[(v * m + u)] = jitter_by(light->u.area.jitter);

                    //canonical[2 * (v * m + u)] = u/(double)m + (v + jitter_by(light->u.area.jitter)) / (n*m);
                    //canonical[2 * (v * m + u) + 1] = v/(double)n + (u + jitter_by(light->u.area.jitter)) / (n*m);
                    //printf("%f %f\n", canonicalx[(v * m + u)], canonicaly[(v * m + u)]);
                }
                //printf("\n");
            }
            //printf("\n");

            // shuffle canonical for x
            for (v = 0; v < n; ++v) {
                k = v + jitter_by(light->u.area.jitter) * (double)(n - v);
                for (u = 0; u < m; ++u) {
                    double_swap(canonicalx + (v * m + u), canonicalx + (k * m + u));
                }
            }

            // shuffle canonical for y
            for (u = 0; u < m; ++u) {
                k = u + jitter_by(light->u.area.jitter) * (double)(m - u);
                for (v = 0; v < n; ++v) {
                    double_swap(canonicaly + (v * m + u), canonicaly + (v * m + k));
                }
            }

            for (v = 0; v < n; ++v) {
                for (u = 0; u < m; ++u) {
                    //printf("%f %f\n", canonicalx[(v * m + u)], canonicaly[(v * m + u)]);
                    area_light_point_on_light(light, u, v, m * canonicalx[(v * m + u)], n * canonicaly[(v * m + u)], point);
                    point_copy(*(itr->points + v * m + u), point);
                }
                //printf("\n");
            }
            //printf("\n");
        }
        light->surface_points_cache = pts;
        light->surface_points_cache_len = AREA_LIGHT_CACHE_SIZE;
        free(canonicalx);
        free(canonicaly);
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
                 Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    area_light(corner, full_uvec, usteps, full_vvec, vsteps, jitter, intensity, l);
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
