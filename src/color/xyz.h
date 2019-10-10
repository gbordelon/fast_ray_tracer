#ifndef XYZ_COLOR
#define XYZ_COLOR

/*
    https://en.wikipedia.org/wiki/SRGB
*/
static double srgb_to_xyz_transform[] = {
    0.412453, 0.357580, 0.1805,
    0.2126, 0.715160, 0.072169,
    0.019334, 0.119193, 0.9505
};

/*
    https://www.cs.rit.edu/~ncs/color/t_convert.html

   [ R ]   [  3.240479 -1.537150 -0.498535 ]   [ X ]
   [ G ] = [ -0.969256  1.875992  0.041556 ] * [ Y ]
   [ B ]   [  0.055648 -0.204043  1.057311 ]   [ Z ]

The range for valid R, G, B values is [0,1]. Note, this matrix has negative coefficients. Some XYZ color may be transformed to RGB values that are negative or greater than one. This means that not all visible colors can be produced using the RGB system.

The inverse transformation matrix is as follows:

   [ X ]   [  0.412453  0.357580  0.180423 ]   [ R ]
   [ Y ] = [  0.212671  0.715160  0.072169 ] * [ G ]
   [ Z ]   [  0.019334  0.119193  0.950227 ]   [ B ]
*/
static double rgb_to_xyz_transform[] = {
    0.412453, 0.357580, 0.180423,
    0.212671, 0.715160, 0.072169,
    0.019334, 0.119193, 0.950227
};

static double xyz_to_rgb_transform[] = {
    3.240479, -1.537150, -0.498535,
    -0.969256, 1.875992, 0.041556,
    0.055648, -0.204043, 1.057311
};

void xyz_to_rgb(const Color xyz, Color rgb);
void xyz_to_srgb(const Color xyz, Color srgb);
void xyz_to_lab(const Color xyz, Color lab);

#endif
