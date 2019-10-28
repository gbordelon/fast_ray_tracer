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
    AREA_LIGHT,
    CIRCLE_LIGHT,
    HEMISPHERE_LIGHT,
    POINT_LIGHT,
    SPOT_LIGHT
};

struct spot_light_fields {
    Point position;
    Vector normal;
    double outer_angle;
    double inner_angle;
};

struct hemisphere_light_fields {
    Point position;
    Vector normal;
};

struct point_light_fields {
    Point position;
};

struct circle_light_fields {
    Point origin;
    Vector normal;
    double radius;
    size_t usteps;
    size_t vsteps;
    bool jitter;
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
        struct area_light_fields area;
        struct circle_light_fields circle;
        struct hemisphere_light_fields hemi;
        struct point_light_fields point;
        struct spot_light_fields spot;
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
void spot_light(Point p, Vector normal, double outer_angle, double inner_angle, Color intensity, Light l);

void
area_light(Point corner,
           Vector full_uvec,
           size_t usteps,
           Vector full_vvec,
           size_t vsteps,
           bool jitter,
           size_t cache_size,
           Color intensity,
           Light l);

void
circle_light(Point origin,
             Vector normal,
             double radius,
             size_t usteps,
             size_t vsteps,
             bool jitter,
             size_t cache_size,
             Color intensity,
             Light l);

#endif
