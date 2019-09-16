#ifndef RENDERER
#define RENDERER

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

typedef struct world {
    Light lights;
    Shape shapes;
    size_t lights_num;
    size_t shapes_num;
} *World;


#endif
