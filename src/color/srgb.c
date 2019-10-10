#include <math.h>

#include "rgb.h"
#include "srgb.h"

void
srgb_to_xyz(const Color srgb, Color xyz)
{
    Color rgb;

    srgb_to_rgb(srgb, rgb);
    rgb_to_xyz(rgb, xyz);
}

void
srgb_to_rgb(const Color srgb, Color rgb)
{
    int i;
    for (i = 0; i < 3; ++i) {
        rgb[i] = srgb[i] <= 0.04045
               ? srgb[i] / 12.92
               : pow((srgb[i] + 0.055) / 1.055, 2.4);
    }
}
