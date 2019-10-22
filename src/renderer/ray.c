#include <stdio.h>

#include "../libs/linalg/linalg.h"
#include "ray.h"

void
ray_array(Point origin, Vector direction, Ray ray)
{
    point_copy(ray->origin, origin);
    vector_copy(ray->direction, direction);
}

void
ray_transform(Ray original, Matrix m, Ray res)
{
    matrix_point_multiply(m, original->origin, res->origin);
    matrix_vector_multiply(m, original->direction, res->direction);
}

int
ray_to_string(char *buf, size_t n, Ray r)
{
    return snprintf(buf, n, "Point: [%f %f %f] Vector: [%f %f %f]",
                    r->origin[0], r->origin[1], r->origin[2],
                    r->direction[0], r->direction[1], r->direction[2]);
}

void
ray_position(Ray ray, double t, Point position)
{
    point_copy(position, ray->origin);
    position[0] += ray->direction[0] * t;
    position[1] += ray->direction[1] * t;
    position[2] += ray->direction[2] * t;
}
