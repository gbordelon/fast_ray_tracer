#ifndef LIGHT
#define LIGHT

#include <stdbool.h>
#include <stdlib.h>

#include "../libs/linalg/linalg.h"
#include "../color/color.h"
#include "../renderer/ray.h"
#include "../renderer/world.h"

struct world;
typedef struct world *World;

enum light_enum {
    POINT_LIGHT,
    AREA_LIGHT,
    HEMISPHERE_LIGHT,
    SPOT_LIGHT
};

struct hemisphere_light_fields {
    Point position;
    Vector normal;
};

struct point_light_fields {
    Point position;
};

struct area_light_fields {
    Point corner;
    Vector uvec;
    size_t usteps;
    Vector vvec;
    size_t vsteps;
    bool jitter;
};

typedef struct light {
    enum light_enum type;
    double intensity[3];
    size_t num_samples;
    size_t num_photons;

    union {
        struct point_light_fields point;
        struct area_light_fields area;
        struct hemisphere_light_fields hemi;
    } u;
    Points surface_points_cache;
    size_t surface_points_cache_len;

    Points (*light_surface_points)(struct light *light);
    double (*intensity_at)(struct light *light, World world, Point point);
    void (*emit_photon)(struct light *light, Ray res);
} *Light;

Light array_of_lights(size_t num);
void point_light(Point p, Color intensity, Light l);
void hemisphere_light(Point p, Vector normal, Color intensity, Light l);
void
area_light(Point corner,
           Vector full_uvec/*vector*/,
           size_t usteps,
           Vector full_vvec/*vector*/,
           size_t vsteps,
           bool jitter,
           Color intensity,
           Light l);

#endif
