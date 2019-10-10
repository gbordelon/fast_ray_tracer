#ifndef WORLD
#define WORLD

#include <stdlib.h>

#include "../libs/photon_map/pm.h"
#include "../light/light.h"
#include "../shapes/shapes.h"
#include "../intersection/intersection.h"

#include "config.h"
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
    PhotonMap *photon_maps;
    Global_config global_config;
} *World;

World world();
void world_copy(const World w, World res);
World default_world();

Intersections intersect_world(const World w, const Ray r);

#endif
