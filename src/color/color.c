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
color_copy(Color to, const Color from)
{
    memcpy(to, from, sizeof(Color));
}
