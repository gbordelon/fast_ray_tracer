#ifndef LIGHT
#define LIGHT

#include <stdbool.h>
#include <stdlib.h>

#include "../libs/linalg/linalg.h"
#include "../color/color.h"
#include "../renderer/world.h"

struct world;
typedef struct world *World;

enum light_enum {
    POINT_LIGHT,
    AREA_LIGHT
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
    double photon_count;

    union {
        struct point_light_fields point;
        struct area_light_fields area;
    } u;
    Points surface_points_cache;
    size_t surface_points_cache_len;

    Points (*light_surface_points)(struct light *light);
    double (*intensity_at)(struct light *light, World world, Point point);
} *Light;

Light array_of_lights(size_t num);
void point_light(Point p, Color intensity, double photon_count, Light l);
void
area_light(Point corner,
           Vector full_uvec/*vector*/,
           size_t usteps,
           Vector full_vvec/*vector*/,
           size_t vsteps,
           bool jitter,
           Color intensity,
           double photon_count,
           Light l);

#endif
