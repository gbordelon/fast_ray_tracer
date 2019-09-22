#ifndef BOUNDING_BOX
#define BOUNDING_BOX

#include <stdbool.h>

#include "shapes.h"

typedef struct ray *Ray;

typedef struct bbox {
    double min[4]; // Point
    double max[4]; // Point
} *Bounding_box;

void bounding_box(Bounding_box res);
Bounding_box bounding_box_alloc();

void bounding_box_free(Bounding_box box);

void bounding_box_add_point(Bounding_box box, Point p);
void bounding_box_add_array(Bounding_box box, double point[4]);
void bounding_box_add_box(Bounding_box box, Bounding_box other);

bool bounding_box_contains_point(Bounding_box box, Point p);
bool bounding_box_contains_array(Bounding_box box, double point[4]);
bool bounding_box_contains_box(Bounding_box box, Bounding_box other);

// this alters box
void bounding_box_transform(Bounding_box box, Matrix m);

bool bounding_box_intersects(Bounding_box, Ray r);

Bounding_box bounding_box_split_bounds(Bounding_box box);

#endif
