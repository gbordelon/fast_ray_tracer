#include "color.h"

void
rgb_to_hsl(Color rgb, Color hsl)
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
hsl_to_rgb(Color hsl, Color rgb)
{
}

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
void
rgb_to_xyz(Color rgb, Color xyz)
{
    int i;
    static double xform[] = {
        0.412453, 0.357580, 0.180423,
        0.212671, 0.715160, 0.072169,
        0.019334, 0.119193, 0.950227
    };

    for (i = 0; i < 3; i++) {
        xyz[i] = xform[i*3+0] * rgb[0] +
                 xform[i*3+1] * rgb[1] +
                 xform[i*3+2] * rgb[2];
    }
}

void
xyz_to_rgb(Color xyz, Color rgb)
{
    int i;
    static double xform[] = {
        3.240479, -1.537150, -0.498535,
        -0.969256, 1.875992, 0.041556,
        0.055648, -0.204043, 1.057311
    };

    for (i = 0; i < 3; i++) {
        rgb[i] = xform[i*3+0] * xyz[0] +
                 xform[i*3+1] * xyz[1] +
                 xform[i*3+2] * xyz[2];
    }
}

/*
    https://en.wikipedia.org/wiki/Illuminant_D65
    https://www.mathworks.com/help/images/ref/rgb2lab.html
*/
static const double tristimulus_2deg[] = {
    0.95047,
    1.00000,
    1.08883
};

static const double tristimulus_10deg[] = {
    0.948110,
    1.000000,
    1.073040
};

static const double *tristimulus = tristimulus_10deg;

void
xyz_to_lab(Color xyz, Color lab)
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

void
lab_to_xyz(Color lab, Color xyz)
{
    double p = (lab[0] + 16.0) / 116.0;
    xyz[0] = tristimulus[0] * pow(p + lab[1] / 500.0, 3.0);
    xyz[1] = tristimulus[1] * pow(p, 3.0);
    xyz[2] = tristimulus[2] * pow(p - lab[2] / 200.0, 3.0);
}

void
rgb_to_lab(Color rgb, Color lab)
{
    Color xyz;
    rgb_to_xyz(rgb, xyz);
    xyz_to_lab(xyz, lab);
}

void
lab_to_rgb(Color lab, Color rgb)
{
    Color xyz;
    lab_to_xyz(lab, xyz);
    xyz_to_rgb(xyz, rgb);
}

int
lab_compare_l(Color l, Color r)
{
    if (l[0] - r[0] < 0) {
        return -1;
    } else if (l[0] - r[0] > 0) {
        return 1;
    }
    return 0;
}

int
lab_compare_a(Color l, Color r)
{
    if (l[1] - r[1] < 0) {
        return -1;
    } else if (l[1] - r[1] > 0) {
        return 1;
    }
    return 0;
}

int
lab_compare_b(Color l, Color r)
{
    if (l[2] - r[2] < 0) {
        return -1;
    } else if (l[2] - r[2] > 0) {
        return 1;
    }
    return 0;
}

void
color_accumulate(Color acc, Color other)
{
    acc[0] += other[0];
    acc[1] += other[1];
    acc[2] += other[2];
}

void
color_scale(Color acc, double scalar)
{
    acc[0] *= scalar;
    acc[1] *= scalar;
    acc[2] *= scalar;
}

int
color_to_string(char *buf, size_t n, Color c)
{
    return snprintf(buf, n, "Color: [%f %f %f]", c[0], c[1], c[2]);
}

void
print_color(Color c)
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



