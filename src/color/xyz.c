#include <math.h>

#include "color.h"
#include "lab.h"
#include "rgb.h"
#include "xyz.h"

void
xyz_to_rgb(const Color xyz, Color rgb)
{
    int i;
    double *xform = xyz_to_rgb_transform;

    for (i = 0; i < 3; i++) {
        rgb[i] = xform[i*3+0] * xyz[0] +
                 xform[i*3+1] * xyz[1] +
                 xform[i*3+2] * xyz[2];
    }
}

void
xyz_to_srgb(const Color xyz, Color srgb)
{
    Color rgb;

    xyz_to_rgb(xyz, rgb);
    rgb_to_srgb(rgb, srgb);
}

void
xyz_to_lab(const Color xyz, Color lab)
{
    double _x = xyz[0] / tristimulus[0];
    double _y = xyz[1] / tristimulus[1];
    double _z = xyz[2] / tristimulus[2];

    double fx = _x > 0.008856
              ? pow(_x, 1.0 / 3.0)
              : 7.787 * _x + 16.0 / 116.0;

    double fy = _y > 0.008856
              ? pow(_y, 1.0 / 3.0)
              : 7.787 * _y + 16.0 / 116.0;

    double fz = _z > 0.008856
              ? pow(_z, 1.0 / 3.0)
              : 7.787 * _z + 16.0 / 116.0;

    lab[0] = _y > 0.008856
           ? 116.0 * pow(_y, 1.0 / 3.0) - 16.0
           : 903.3 * _y;

    lab[1] = 500.0 * (fx - fy);
    lab[2] = 200.0 * (fy - fz);
}
