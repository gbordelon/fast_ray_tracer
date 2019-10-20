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
#include "ray.h"

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
    struct ray photon_ray;
    Color photon_power;
    Color over_Ka;
    Color over_Kd;
    Color over_Ks;
    Color over_refl;
} *Computations;

struct container;

void shade_hit(World w, Computations comps, size_t remaining, Color res);
bool is_shadowed(World w, Point light_position, Point pt);

double schlick(Computations comps);
void prepare_computations(Intersection i, Ray r, Color photon_power, Intersections xs, Computations res, struct container *container);
Canvas render(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter);
Canvas render_multi(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter);

#endif
