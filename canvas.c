#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>

#include "canvas.h"
#include "linalg.h"

#define header_len 32

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
   [ B ]   [  0.055648 -0.204043  1.057311 ]   [ Z ].

The range for valid R, G, B values is [0,1]. Note, this matrix has negative coefficients. Some XYZ color may be transformed to RGB values that are negative or greater than one. This means that not all visible colors can be produced using the RGB system.

The inverse transformation matrix is as follows:

   [ X ]   [  0.412453  0.357580  0.180423 ]   [ R ] **
   [ Y ] = [  0.212671  0.715160  0.072169 ] * [ G ]
   [ Z ]   [  0.019334  0.119193  0.950227 ]   [ B ].
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
color_copy(Color to, Color from)
{
    memcpy(to, from, sizeof(Color));
}


Canvas
canvas_alloc(size_t width, size_t height)
{
    Canvas c = (Canvas) malloc(sizeof(struct canvas));
    // null check c
    c->arr = (Color *) malloc(width * height * sizeof(Color));
    // null check c->arr

    c->width = width;
    c->height = height;
    c->depth = 3;

    pthread_mutex_init(&(c->write_lock), NULL);

    return c;
}

/*
 * len should be width * height * 3 + 1
 * 3 because one byte for each of rgb
 * 1 because trailing newline
 */
Ppm
ppm_alloc(size_t len)
{
    Ppm p = (Ppm) malloc(sizeof(struct ppm_struct));
    // null check p
    p->arr = (unsigned char*) malloc(len * sizeof(unsigned char));
    // null check arr
    p->len = len;

    return p;
}

void
canvas_free(Canvas c)
{
    if (c != NULL) {
        if (c->arr != NULL) {
            free(c->arr);
        }
        free(c);
    }
}

void
ppm_free(Ppm p)
{
    if (p != NULL) {
        if (p->arr != NULL) {
            free(p->arr);
        }
        free(p);
    }
}

/*
 * row and col do not take into account 3 doubles for each entry
 * e.g. 4x4 canvas
 * 00 -> 0
 * 01 -> 3
 * 02 -> 6
 * 03 -> 9
 *
 * 10 -> 12
 * 11 -> 15
 * 12 -> 18
 * 13 -> 21
 *
 * 20 -> 24
 * 21 -> 27
 * 22 -> 30
 * 23 -> 33
 *
 * 30 -> 36
 * 31 -> 39
 * 32 -> 42
 * 33 -> 45
 */
void
canvas_write_pixel(Canvas c, int col, int row, Color color)
{
    // null check c
    pthread_mutex_lock(&c->write_lock);
    color_copy(*(c->arr + row * c->width + col), color);
    pthread_mutex_unlock(&c->write_lock);
}

void
canvas_pixel_at(Canvas c, int col, int row, Color res)
{
    // null check c
    color_copy(res, *(c->arr + row * c->width + col));
}

Ppm
construct_ppm(Canvas c, bool use_scaling)
{
    int n, i, out_len;
    unsigned char *out, *buf;
    uint16_t r_scaled, g_scaled, b_scaled;
    double r_inverse, g_inverse, b_inverse;
    Color lab_max, lab_tmp, rgb_max, rgb_tmp;
    Color *cur_val;
    Ppm ppm;

    buf = (unsigned char *) malloc(header_len * sizeof(unsigned char));
    if (buf == NULL) {
        // error
        return NULL;
    }

    n = snprintf((char *)buf, header_len, "P6\n%zu %zu\n65535\n", c->width, c->height);
    // error check n

    out_len = c->width * c->height * c->depth * 2 + n + 1;
    ppm = ppm_alloc(out_len);

    if (ppm == NULL) {
        free(buf);
        // error
    }
    out = ppm->arr;

    memcpy(out, buf, n * sizeof(unsigned char));
    out += n;

    free(buf);

    rgb_max[0] = 1.0;
    rgb_max[1] = 1.0;
    rgb_max[2] = 1.0;

    // iterate over all colors in c->arr and find the max
    if (use_scaling) {
        // iterate to find max rgb for scaling to provide to max_lab
        for (i = 0; i < c->width * c->height; i++) {
            if (*(c->arr + i)[0] > rgb_max[0]) {
                rgb_max[0] = *(c->arr + i)[0];
            }
            if (*(c->arr + i)[1] > rgb_max[1]) {
                rgb_max[1] = *(c->arr + i)[1];
            }
            if (*(c->arr + i)[2] > rgb_max[2]) {
                rgb_max[2] = *(c->arr + i)[2];
            }
        }
        lab_max[0] = 0;
        lab_max[1] = 0;
        lab_max[2] = 0;

        // iterate to find max L*
        for (i = 0; i < c->width * c->height; i++) {
            color_copy(rgb_tmp, *(c->arr + i));
            rgb_tmp[0] /= rgb_max[0];
            rgb_tmp[1] /= rgb_max[1];
            rgb_tmp[2] /= rgb_max[2];

            rgb_to_lab(rgb_tmp, lab_tmp);
            if (lab_compare_l(lab_tmp, lab_max) > 0) {
                lab_max[0] = lab_tmp[0];
            }
        }

        printf("lab_max: %f %f %f\n", lab_max[0], lab_max[1], lab_max[2]);
        lab_to_rgb(lab_max, rgb_max);
    } else {
        rgb_max[0] = 1.0;
        rgb_max[1] = 1.0;
        rgb_max[2] = 1.0;
    }

    printf("rgb_max: %f %f %f\n", rgb_max[0], rgb_max[1], rgb_max[2]); 
    r_inverse = 65535.0 / rgb_max[0];
    g_inverse = 65535.0 / rgb_max[1];
    b_inverse = 65535.0 / rgb_max[2];

    for (cur_val = c->arr; n < out_len - 1; n += 6, cur_val++) {
        if ((*cur_val)[0] >= rgb_max[0]) {
            r_scaled = 65535;
        } else if ((*cur_val)[0] <= 0) {
            r_scaled = 0;
        } else {
            r_scaled = (uint16_t)floor((*cur_val)[0] * r_inverse);
        }
        *out++ = (r_scaled & 0xFF00) >> 8;
        *out++ = (r_scaled & 0x00FF);

        if ((*cur_val)[1] >= rgb_max[1]) {
            g_scaled = 65535;
        } else if ((*cur_val)[1] <= 0) {
            g_scaled = 0;
        } else {
            g_scaled = (uint16_t)floor((*cur_val)[1] * g_inverse);
        }
        *out++ = (g_scaled & 0xFF00) >> 8;
        *out++ = (g_scaled & 0x00FF);

        if ((*cur_val)[2] >= rgb_max[2]) {
            b_scaled = 65535;
        } else if ((*cur_val)[2] <= 0) {
            b_scaled = 0;
        } else {
            b_scaled = (uint16_t)floor((*cur_val)[2] * b_inverse);
        }
        *out++ = (b_scaled & 0xFF00) >> 8;
        *out++ = (b_scaled & 0x00FF);
    }
    *out = '\n';

    return ppm;
}

Canvas
construct_canvas_from_ppm_file(const char * file_path)
{
    char buf[32];
    size_t width, height, max_val;
    int i, total_rgb_count;
    unsigned int tmp;
    Color *out;

    FILE *file = fopen(file_path, "r");
    if (file == NULL) {
        printf("Error opening file %s", file_path);
        return NULL;
    }
    fscanf(file, "%s", buf); // P6
    fscanf(file, "%zu", &width); // width
    fscanf(file, "%zu", &height); // height
    fscanf(file, "%zu", &max_val); // max val

    Canvas c = canvas_alloc(width, height);
    if (c == NULL) {
        fclose(file);
        return NULL;
    }

    total_rgb_count = width * height * 3;
    for (i = 0, out = c->arr; i < total_rgb_count; i++, out++) {
        fscanf(file, "%u", &tmp);
        *out[0] = ((double) tmp) / (double)max_val;
        fscanf(file, "%u", &tmp);
        *out[1] = ((double) tmp) / (double)max_val;
        fscanf(file, "%u", &tmp);
        *out[2] = ((double) tmp) / (double)max_val;
    }

    fclose(file);

    return c;
}
