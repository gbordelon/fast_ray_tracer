#include <stdlib.h>
#include <math.h>

#include "bounding_box.h"

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
bounding_box_free(Bounding_box box)
{
    if (box != NULL) {
        free(box);
    }
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

    box->min[0] = min_x;
    box->min[1] = min_y;
    box->min[2] = min_z;

    box->max[0] = max_x;
    box->max[1] = max_y;
    box->max[2] = max_z;
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
    return bounding_box_contains_array(box, p->arr);
}

bool
bounding_box_contains_box(Bounding_box box, Bounding_box other)
{
    return (box == other)
        || (bounding_box_contains_array(box, other->min)
           && bounding_box_contains_array(box, other->max));
}

Bounding_box
bounding_box_transform(Bounding_box box, Matrix m)
{
    double points[8][4] = {
        {box->min[0], box->min[1], box->min[2], 1.0 },
        {box->min[0], box->min[1], box->max[2], 1.0 },
        {box->min[0], box->max[1], box->min[2], 1.0 },
        {box->min[0], box->max[1], box->max[2], 1.0 },
        {box->max[0], box->min[1], box->min[2], 1.0 },
        {box->max[0], box->min[1], box->max[2], 1.0 },
        {box->max[0], box->max[1], box->min[2], 1.0 },
        {box->max[0], box->max[1], box->max[2], 1.0 },
    };

    Bounding_box tr = bounding_box_alloc();
    double arr[4];
    int i;
    for (i = 0; i < 8; i++) {
        matrix_array_multiply(m, points[i], arr);
        bounding_box_add_array(tr, arr);
    }

    return tr;
}

struct two_doubles {
    double a[2];
};

struct two_doubles
bbox_check_axis(double origin, double direction, double min, double max)
{
    struct two_doubles min_max;

    double tmin_numerator = min - origin;
    double tmax_numerator = max - origin;
    double tmin, tmax;

    if (fabs(direction) >= EPSILON) {
        tmin = tmin_numerator / direction;
        tmax = tmax_numerator / direction;
    } else {
        tmin = tmin_numerator * INFINITY;
        if (isnan(tmin)) {
            tmin = INFINITY;
            if (tmin_numerator < 0) {
                tmin = -INFINITY;
            }
        }
        tmax = tmax_numerator * INFINITY;
        if (isnan(tmax)) {
            tmax = INFINITY;
            if (tmax_numerator < 0) {
                tmax = -INFINITY;
            }
        }
    }

    if (tmin > tmax) {
        min_max.a[0] = tmax;
        min_max.a[1] = tmin;
    } else {
        min_max.a[0] = tmin;
        min_max.a[1] = tmax;
    }

    return min_max;

}

bool
bounding_box_intersects(Bounding_box box, Ray r)
{
    struct two_doubles xtmin_xtmax = bbox_check_axis(r->origin[0], r->direction[0], box->min[0], box->max[0]);
    struct two_doubles ytmin_ytmax = bbox_check_axis(r->origin[1], r->direction[1], box->min[1], box->max[1]);
    struct two_doubles ztmin_ztmax = bbox_check_axis(r->origin[2], r->direction[2], box->min[2], box->max[2]);

    double tmin = fmax(fmax(xtmin_xtmax.a[0], ytmin_ytmax.a[0]), ztmin_ztmax.a[0]);
    double tmax = fmin(fmin(xtmin_xtmax.a[1], ytmin_ytmax.a[1]), ztmin_ztmax.a[1]);

    return tmin <= tmax;
}

Bounding_box
bounding_box_split_bounds(Bounding_box box)
{
    double dx = fabs(box->max[0] - box->min[0]);
    double dy = fabs(box->max[1] - box->min[1]);
    double dz = fabs(box->max[2] - box->min[2]);
    double greatest = fmax(fmax(dx,dy),dz);

    double x0 = box->min[0];
    double y0 = box->min[1];
    double z0 = box->min[2];

    double x1 = box->max[0];
    double y1 = box->max[1];
    double z1 = box->max[2];

    if (equal(greatest, dx)) {
        x0 = x1 = x0 + dx / 2.0;
    } else if (equal(greatest, dy)) {
        y0 = y1 = y0 + dy / 2.0;
    } else {
        z0 = z1 = z0 + dz / 2.0;
    }

    Bounding_box boxes = (Bounding_box) malloc(2 * sizeof(struct bbox));
    boxes->min[0] = box->min[0];
    boxes->min[1] = box->min[1];
    boxes->min[2] = box->min[2];
    boxes->max[0] = x1;
    boxes->max[1] = y1;
    boxes->max[2] = z1;

    (boxes+1)->min[0] = x0;
    (boxes+1)->min[1] = y0;
    (boxes+1)->min[2] = z0;
    (boxes+1)->max[0] = box->max[0];
    (boxes+1)->max[1] = box->max[1];
    (boxes+1)->max[2] = box->max[2];

    return boxes;
}

int
bounding_box_to_string(char * buf, size_t n, Bounding_box box)
{
    return snprintf(buf, n, "Box: [%f %f %f] [%f %f %f]",
                    box->min[0], box->min[1], box->min[2],
                    box->max[0], box->max[1], box->max[2]);
}

