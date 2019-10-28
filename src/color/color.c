#include <stdio.h>
#include <string.h>

#include "color.h"

void
color_accumulate(Color acc, const Color other)
{
    acc[0] += other[0];
    acc[1] += other[1];
    acc[2] += other[2];
}

void
color_scale(Color acc, const double scalar)
{
    acc[0] *= scalar;
    acc[1] *= scalar;
    acc[2] *= scalar;
}

int
color_to_string(char *buf, size_t n, const Color c)
{
    return snprintf(buf, n, "Color: [%f %f %f]", c[0], c[1], c[2]);
}

void
print_color(const Color c)
{
    static char buf[256];
    color_to_string(buf, 256, c);
    printf("%s\n",buf);
}

void
print_color_triple(const ColorTriple c)
{
    print_color(ambient_from_triple(c));
    print_color(diffuse_from_triple(c));
    print_color(specular_from_triple(c));
    printf("\n");
}

void
color_copy(Color to, const Color from)
{
    memcpy(to, from, sizeof(Color));
}

void
color_triple_copy(ColorTriple to, const ColorTriple from)
{
    memcpy(to, from, sizeof(ColorTriple));
}

void
color_average(Color c1, Color c2, Color res)
{
    res[0] = (c1[0] + c2[0]) / 2.0;
    res[1] = (c1[1] + c2[1]) / 2.0;
    res[2] = (c1[2] + c2[2]) / 2.0;
}

void
color_triple_average(Color c1, Color c2, Color res)
{
    color_average(ambient_from_triple(c1), ambient_from_triple(c2), ambient_from_triple(res));
    color_average(diffuse_from_triple(c1), diffuse_from_triple(c2), diffuse_from_triple(res));
    color_average(specular_from_triple(c1), specular_from_triple(c2), specular_from_triple(res));
}
