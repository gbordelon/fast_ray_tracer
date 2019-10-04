#ifndef WORLD
#define WORLD

#include <stdlib.h>

#include "../light/light.h"
#include "../shapes/shapes.h"
#include "../intersection/intersection.h"

#include "renderer.h"

struct light;
struct ray;

typedef struct light *Light;
typedef struct ray *Ray;

typedef struct world {
    Light lights;
    Shape shapes;
    Intersections xs;
    size_t lights_num;
    size_t shapes_num;
} *World;

World world();
void world_copy(const World w, World res);
World default_world();

Intersections intersect_world(const World w, const Ray r);

#endif
