#ifndef SRGB_COLOR
#define SRGB_COLOR

#include "color.h"

void srgb_to_xyz(const Color rgb, Color xyz);
void srgb_to_rgb(const Color srgb, Color rgb);

#endif

