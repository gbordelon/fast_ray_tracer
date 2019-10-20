#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <stdbool.h>

#include "../libs/linalg/linalg.h"

typedef struct ray *Ray;

typedef struct bbox {
    Point min;
    Point max;
} Bounding_box;

#include "shapes.h"

void bounding_box(Bounding_box *res);

void bounding_box_add_point(Bounding_box *box, Point p);
void bounding_box_add_array(Bounding_box *box, double point[4]);
void bounding_box_add_box(Bounding_box *box, Bounding_box *other);

bool bounding_box_contains_point(Bounding_box *box, Point p);
bool bounding_box_contains_array(Bounding_box *box, double point[4]);
bool bounding_box_contains_box(Bounding_box *box, Bounding_box *other);

void bounding_box_transform(Bounding_box *box, Matrix m, Bounding_box *res);

bool bounding_box_intersects(Bounding_box *box, Ray r);

void bounding_box_split_bounds(Bounding_box *box, Bounding_box *left_res, Bounding_box *right_res);

int bounding_box_to_string(char * buf, size_t n, Bounding_box *box);

#endif
