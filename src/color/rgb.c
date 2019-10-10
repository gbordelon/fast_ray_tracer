#include <math.h>

#include "../libs/linalg/linalg.h"

#include "color.h"
#include "rgb.h"
#include "xyz.h"

void
rgb_to_rgb(const Color from, Color to)
{
    color_copy(to, from);
}

void
rgb_to_hsl(const Color rgb, Color hsl)
{
    double r = rgb[0];
    double g = rgb[1];
    double b = rgb[2];
    double max = fmax(fmax(r, g), b);
    double min = fmin(fmin(r, g), b);

    hsl[2] = (max + min) / 2.0;

    if (hsl[2] < 0.5) {
        hsl[1] = (max - min) / (max + min);
    } else {
        hsl[1] = (max - min) / (2.0 - max - min);
    }

    if (equal(max, r)) {
        hsl[0] = (g - b) / (max - min);
    } else if (equal(max, g)) {
        hsl[0] = 2.0 + (b - r) / (max - min);
    } else {
        hsl[0] = 4.0 + (r - g) / (max - min);
    }
    hsl[0] *= 60;
    if (hsl[0] < 0) {
        hsl[0] += 360.0;
    }
}

void
rgb_to_xyz(const Color rgb, Color xyz)
{
    int i;
    double *xform = rgb_to_xyz_transform;

    for (i = 0; i < 3; i++) {
        xyz[i] = xform[i*3+0] * rgb[0] +
                 xform[i*3+1] * rgb[1] +
                 xform[i*3+2] * rgb[2];
    }
}

void
rgb_to_lab(const Color rgb, Color lab)
{
    Color xyz;
    rgb_to_xyz(rgb, xyz);
    xyz_to_lab(xyz, lab);
}

void
rgb_to_srgb(const Color rgb, Color srgb)
{
    int i;

    // gamma correct the rgb to srgb
    for (i = 0; i < 3; ++i) {
        srgb[i] = rgb[i] < 0.0031308
                ? rgb[i] * 12.92 
                : (1.055 * pow(rgb[i], 1.0 / 2.4) - 0.055);
    }
}
