#ifndef RENDERER
#define RENDERER

#include <stdbool.h>
#include <stdlib.h>

#include "../libs/linalg/linalg.h"
#include "../color/color.h"
#include "../light/light.h"
#include "../shapes/shapes.h"
#include "../intersection/intersection.h"

#include "world.h"
#include "camera.h"

typedef struct ray {
    Point origin;
    Vector direction;
} *Ray;

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

struct container {
    Shape *shapes;
    size_t size;
};

void ray_array(Point origin, Vector direction, Ray ray);
int ray_to_string(char *s, size_t n, Ray r);
void ray_transform(Ray original, Matrix m, Ray res);

void shade_hit(World w, Computations comps, size_t remaining, Color res, struct container *container);
bool is_shadowed(World w, Point light_position, Point pt);

void prepare_computations(Intersection i, Ray r, Intersections xs, Computations res, struct container *container);
Canvas render(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter);
Canvas render_multi(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter, size_t num_threads);

#endif
