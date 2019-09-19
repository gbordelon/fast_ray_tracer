#ifndef RENDERER
#define RENDERER

#include <stdbool.h>
#include <stdlib.h>

#include "canvas.h"
#include "linalg.h"
#include "shapes.h"

// World
// Computations
// Camera
// Light
// PointLight
// AreaLight
// prepare computations
// intersect world
// schlick
// refracted color
// relfected color
// reflect
// shade hit
// is shadowed
// color at
// view transform
// ray for pixel
// render
// normal at
// lighting

typedef struct world *World;

enum light_enum {
    POINT_LIGHT,
    AREA_LIGHT
};

struct point_light_fields {
    double position[4];
};

struct area_light_fields {
    double corner[4];
    double uvec[4];
    size_t usteps;
    double vvec[4];
    size_t vsteps;
    bool jitter;
};

typedef struct light {
    enum light_enum type;
    double intensity[3];
    size_t num_samples;
    union {
        struct point_light_fields point;
        struct area_light_fields area;
    } u;
    Points surface_points_cache;
    size_t surface_points_cache_len;

    Points (*light_surface_points)(struct light *light);
    double (*intensity_at)(struct light *light, World world, Point point);
} *Light;

typedef struct computations {
    double t;
    double n1;
    double n2;
    bool inside;
    Shape obj;
    Point p;
    Point over_point;
    Point under_point;
    Vector eyev;
    Vector normalv;
    Vector reflectv;
} *Computations;

typedef struct world {
    Light lights;
    Shape shapes;
    size_t lights_num;
    size_t shapes_num;
} *World;

typedef struct camera {
    size_t hsize;
    size_t vsize;
    double field_of_view;
    double half_width;
    double half_height;
    double pixel_size;
    Matrix transform;
    Matrix transform_inverse;
} *Camera;

Camera camera(size_t hsize, size_t vsize, double field_of_view, Matrix transform);
Matrix view_transform(Point fr, Point to, Vector up);
void camera_set_transform(Camera c, Matrix m);

World world();
World default_world();


// pixel x and pixel y
//Ray ray_for_pixel(Camera cam, size_t px, size_t py, double x_offset, double y_offset);
Color color_at(World w, Ray r, size_t remaining);
Color shade_hit(World w, Computations comps, size_t remaining);
bool is_shadowed(World w, double light_position[4], Point pt);
Intersections intersect_world(World w, Ray r);
Computations prepare_computations(Intersection i, Ray r, Intersections xs);
Canvas render(Camera cam, World w);

#endif
