#include <stdlib.h>
#include <math.h>

#include "bounding_box.h"

/*
typedef struct bbox {
    double min[4]; // Point
    double max[4]; // Point
} *Bounding_Box;

void bounding_box(Bounding_box res);
Bounding_Box bounding_box_alloc();

void bounding_box_add_point(Bounding_box box, Point p);
void bounding_box_add_array(Bounding_box box, double point[4]);
void bounding_box_add_box(Bounding_box box, Bounding_box other);

bool bounding_box_contains_point(Bounding_box box, Point p);
bool bounding_box_contains_array(Bounding_box box, double point[4]);
bool bounding_box_contains_box(Bounding_box box, Bounding_box other);

void bounding_box_transform(Bounding_box box, Matrix m);

Intersections intersect(Bounding_box, Ray r);

Bounding_box bounding_box_split_bounds(Bounding_box box);
*/

void
bounding_box(Bounding_box box)
{
    box->min[0] = INFINITY;
    box->min[1] = INFINITY;
    box->min[2] = INFINITY;
    box->min[3] = 1.0;

    box->max[0] = -INFINITY;
    box->max[1] = -INFINITY;
    box->max[2] = -INFINITY;
    box->max[3] = 1.0;
}

Bounding_box
bounding_box_alloc()
{
    Bounding_box box = (Bounding_box) malloc(sizeof(struct bbox));
    bounding_box(box);
    return box;
}

void
bounding_box_add_array(Bounding_box box, double point[4])
{
    double min_x = box->min[0];
    double min_y = box->min[1];
    double min_z = box->min[2];
    double max_x = box->max[0];
    double max_y = box->max[1];
    double max_z = box->max[2];

    if (point[0] < box->min[0]) {
        min_x = point[0];
    }
    if (point[1] < box->min[1]) {
        min_y = point[1];
    }
    if (point[2] < box->min[2]) {
        min_z = point[2];
    }
    if (point[0] > box->max[0]) {
        max_x = point[0];
    }
    if (point[1] > box->max[1]) {
        max_y = point[1];
    }
    if (point[2] > box->max[2]) {
        max_z = point[2];
    }

    box->min[0] = min_x
    box->min[1] = min_y
    box->min[2] = min_z

    box->max[0] = max_x
    box->max[1] = max_y
    box->max[2] = max_z
}

void
bounding_box_add_point(Bounding_box box, Point p)
{
    bounding_box_add_array(box, p->arr);
}

void
bounding_box_add_box(Bounding_box box, Bounding_box other)
{
    if (box != other) {
        bounding_box_add_array(box, other->min);
        bounding_box_add_array(box, other->max);
    }
}

bool
bounding_box_contains_array(Bounding_box box, double point[4])
{
    return box->min[0] <= point[0] && point[0] <= box->max[0]
        && box->min[1] <= point[1] && point[1] <= box->max[1]
        && box->min[2] <= point[2] && point[2] <= box->max[2];
}

bool
bounding_box_contains_point(Bounding_box box, Point p)
{
    bounding_box_contains_array(box, p->arr);
}

bool
bounding_box_contains_box(Bounding_box box, Bounding_box other)
{
    return (box == other)
        || (bounding_box_contains_array(box, other->min)
           && bounding_box_contains_array(box, other->max));
}

void
bounding_box_transform(Bounding_box box, Matrix m)
{
}


struct two_doubles {
    double a[2];
};

struct two_doubles
check_axis(double origin, double direction, double min, double max)
{
    struct tow_doubles min_max;
/*
 *  TODO
 */
    return min_max;
}

Intersections
intersect(Bounding_box, Ray r)
{
}

Bounding_box
bounding_box_split_bounds(Bounding_box box)
{
}

