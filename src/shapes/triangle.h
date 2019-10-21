#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "../libs/linalg/linalg.h"

#include "shapes.h"

Shape triangle_point_alloc(Point p1, Point p2, Point p3);
Shape triangle_array_alloc(Point p1, Point p2, Point p3);
Shape smooth_triangle_alloc(Point p1, Point p2, Point p3,
                            Vector n1, Vector n2, Vector n3);

void triangle(Shape s, Point p1, Point p2, Point p3);

void smooth_triangle(Shape s, Point p1, Point p2, Point p3,
                     Vector n1, Vector n2, Vector n3);

void convert_triangle_to_smooth_triangle(Shape s, Vector n1, Vector n2, Vector n3);

#endif
