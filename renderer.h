#ifndef RENDERER
#define RENDERER

#include <stdbool.h>
#include <stdlib.h>

#include "canvas.h"
#include "linalg.h"
#include "shapes.h"

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
    Intersections xs;
} *World;

enum aperture_type {
    CIRCULAR_APERTURE,
    CROSS_APERTURE,
    DIAMOND_APERTURE,
    DOUBLE_CIRCLE_APERTURE,
    HEXAGONAL_APERTURE,
    PENTAGONAL_APERTURE,
    POINT_APERTURE,
    SQUARE_APERTURE,
    OCTAGONAL_APERTURE,
};

struct circle_aperture_args {
    double r1;
};

struct cross_aperture_args {
    double x1;
    double x2;
    double y1;
    double y2;
};

struct diamond_aperture_args {
    double b1;
    double b2;
    double b3;
    double b4;
};

struct double_circle_aperture_args {
    double r1;
    double r2;
};

struct aperture {
    enum aperture_type type;
    double size;
    bool jitter;
    union {
        struct circle_aperture_args circle;
        struct cross_aperture_args cross;
        struct diamond_aperture_args diamond;
        struct double_circle_aperture_args double_circle;
    } u;
};

typedef struct camera {
    size_t hsize;
    size_t vsize;
    size_t sample_num;
    double field_of_view;
    double canvas_distance;
    struct aperture aperture;
    double half_width;
    double half_height;
    double pixel_size;
    Matrix transform;
    Matrix transform_inverse;
} *Camera;

Light array_of_lights(size_t num);
void point_light(Point p, Color intensity, Light l);
void
area_light(Point corner,
           Vector full_uvec/*vector*/,
           size_t usteps,
           Vector full_vvec/*vector*/,
           size_t vsteps,
           bool jitter,
           Color intensity,
           Light l);

Camera camera(size_t hsize, size_t vsize, double field_of_view, double canvas_distance, struct aperture aperture, size_t sample_num, Matrix transform);
void view_transform(Point fr, Point to, Vector up, Matrix res);
void camera_set_transform(Camera c, Matrix m);

World world();
World default_world();


void color_at(World w, Ray r, size_t remaining, Color res);
void shade_hit(World w, Computations comps, size_t remaining, Color res);
bool is_shadowed(World w, Point light_position, Point pt);
Intersections intersect_world(World w, Ray r);

void intersections_sort(Intersections xs);
void intersections_reverse(Intersections xs);

void prepare_computations(Intersection i, Ray r, Intersections xs, Computations res);
Canvas render(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter);
Canvas render_multi(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter);

#endif
