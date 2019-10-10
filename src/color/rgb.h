#ifndef RGB_COLOR
#define RGB_COLOR

#include "color.h"

void rgb_to_rgb(const Color from, Color to);
void rgb_to_hsl(const Color rgb, Color hsl);
void rgb_to_xyz(const Color rgb, Color xyz);
void rgb_to_lab(const Color rgb, Color lab);
void rgb_to_srgb(const Color rgb, Color srgb);


#endif
