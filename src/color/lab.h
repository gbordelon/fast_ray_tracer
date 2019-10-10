#ifndef LAB_COLOR
#define LAB_COLOR

#include "color.h"

/*
    https://en.wikipedia.org/wiki/Illuminant_D65
    https://www.mathworks.com/help/images/ref/rgb2lab.html
*/
#define TWO_DEGREES
#ifdef TWO_DEGREES
// 2 degrees
static const double tristimulus[] = {
    0.95047,
    1.00000,
    1.08883
};
#else
// 10 degrees
static const double tristimulus[] = {
    0.948110,
    1.000000,
    1.073040
};
#endif

void lab_to_xyz(const Color lab, Color xyz);
void lab_to_rgb(const Color lab, Color rgb);

int lab_compare_l(const Color l, const Color r);
int lab_compare_a(const Color l, const Color r);
int lab_compare_b(const Color l, const Color r);

#endif
