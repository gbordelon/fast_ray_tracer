#ifndef TRIANGLE
#define TRIANGLE

#include "shapes.h"

Shape triangle_point_alloc(Point p1, Point p2, Point p3);
Shape triangle_array_alloc(double p1[4], double p2[4], double p3[4]);
void triangle(Shape s, double p1[4], double p2[4], double p3[4]);

void smooth_triangle(Shape s, double p1[4], double p2[4], double p3[4],
                     double n1[4], double n2[4], double n3[4]);
Shape smooth_triangle_alloc(double p1[4], double p2[4], double p3[4],
                            double n1[4], double n2[4], double n3[4]);


#endif
