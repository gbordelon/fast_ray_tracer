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
rgb_to_hsl(double rgb[3], double hsl[3])
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
rgb_to_xyz(double rgb[3], double xyz[3])
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
xyz_to_rgb(double xyz[3], double rgb[3])
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
*/
static const double tristimulus_2deg[] = {
    95.047,
    100.00,
    108.883
};

static const double tristimulus_10deg[] = {
    94.8110,
    100.00,
    107.304
};

static const double *tristimulus = tristimulus_10deg;

void
xyz_to_lab(double xyz[3], double lab[3])
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
lab_to_xyz(double lab[3], double xyz[3])
{
    double p = (lab[0] + 16.0) / 116.0;
    xyz[0] = tristimulus[0] * pow(p + lab[1] / 500.0, 3.0);
    xyz[1] = tristimulus[1] * pow(p, 3.0);
    xyz[2] = tristimulus[2] * pow(p - lab[2] / 200.0, 3.0);
}

void
rgb_to_lab(double rgb[3], double lab[3])
{
    double xyz[3];
    rgb_to_xyz(rgb, xyz);
    xyz_to_lab(xyz, lab);
}

void
lab_to_rgb(double lab[3], double rgb[3])
{
    double xyz[3];
    lab_to_xyz(lab, xyz);
    xyz_to_rgb(xyz, rgb);
}

int
lab_compare_l(double l[3], double r[3])
{
    if (l[0] - r[0] < 0) {
        return -1;
    } else if (l[0] - r[0] > 0) {
        return 1;
    }
    return 0;
}

int
lab_compare_a(double l[3], double r[3])
{
    if (l[1] - r[1] < 0) {
        return -1;
    } else if (l[1] - r[1] > 0) {
        return 1;
    }
    return 0;
}

int
lab_compare_b(double l[3], double r[3])
{
    if (l[2] - r[2] < 0) {
        return -1;
    } else if (l[2] - r[2] > 0) {
        return 1;
    }
    return 0;
}


Color
color_default()
{
    Color c = (Color) malloc(sizeof(struct color));
    // null check c
    return c;
}

Color
color(double r, double g, double b)
{
    Color c = color_default();
    c->arr[0] = r;
    c->arr[1] = g;
    c->arr[2] = b;
    return c;
}

void
color_accumulate(Color acc, Color other)
{
    acc->arr[0] += other->arr[0];
    acc->arr[1] += other->arr[1];
    acc->arr[2] += other->arr[2];
}

void
color_scale(Color acc, double scalar)
{
    acc->arr[0] *= scalar;
    acc->arr[1] *= scalar;
    acc->arr[2] *= scalar;
}

int
color_to_string(char *buf, size_t n, Color c)
{
    return snprintf(buf, n, "Color: [%f %f %f]", c->arr[0], c->arr[1], c->arr[2]);
}

Canvas
canvas_alloc(size_t width, size_t height)
{
    Canvas c = (Canvas) malloc(sizeof(struct canvas));
    // null check c
    c->arr = (double *) malloc(width * height * 3 * sizeof(double));
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

void
color_free(Color c)
{
    if (c != NULL) {
        free(c);
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
    // null check color
    pthread_mutex_lock(&c->write_lock);
    memcpy(c->arr + c->depth * (row * c->width + col), color->arr, sizeof(color->arr));
    pthread_mutex_unlock(&c->write_lock);
}

void
canvas_pixel_at(Canvas c, int col, int row, Color res)
{
    // null check c
    // null check res
    memcpy(res->arr, c->arr + c->depth * (row * c->width + col), sizeof(res->arr));
}

Color
canvas_pixel_at_alloc(Canvas c, int col, int row)
{
    Color color = (Color) malloc(sizeof(struct color));
    // null check color
    canvas_pixel_at(c, row, col, color);
    return color;
}

Ppm
construct_ppm(Canvas c, bool use_scaling)
{
    int n, i, out_len;
    unsigned char *out, *buf;
    uint16_t r_scaled, g_scaled, b_scaled;
    double *cur_val, r_inverse, g_inverse, b_inverse;
    double lab_max[3], lab_tmp[3], rgb_max[3], rgb_tmp[3];
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
        for (i = 0; i < c->width * c->height * c->depth; i += c->depth) {
            if (*(c->arr + i) > rgb_max[0]) {
                rgb_max[0] = *(c->arr + i);
            }
            if (*(c->arr + i + 1) > rgb_max[1]) {
                rgb_max[1] = *(c->arr + i + 1);
            }
            if (*(c->arr + i + 2) > rgb_max[2]) {
                rgb_max[2] = *(c->arr + i + 2);
            }
        }
/*
        lab_max[0] = 0;
        lab_max[1] = 0;
        lab_max[2] = 0;

        // iterate to find max L*
        for (i = 0; i < c->width * c->height * c->depth; i += c->depth) {
            memcpy(rgb_tmp, c->arr + i, 3 * sizeof(double));
            rgb_tmp[0] /= rgb_max[0];
            rgb_tmp[1] /= rgb_max[1];
            rgb_tmp[2] /= rgb_max[2];

            rgb_to_lab(rgb_tmp, lab_tmp);
            if (lab_compare_l(lab_tmp, lab_max) > 0) {
                lab_max[0] = lab_tmp[0];
            }
        }

        printf("lab_max: %f %f %f\n", lab_max[0], lab_max[1], lab_max[2]);
        //lab_to_rgb(lab_max, rgb_max);
*/
    } else {
        rgb_max[0] = 1.0;
        rgb_max[1] = 1.0;
        rgb_max[2] = 1.0;
    }

    printf("rgb_max: %f %f %f\n", rgb_max[0], rgb_max[1], rgb_max[2]); 
    r_inverse = 65535.0 / rgb_max[0];
    g_inverse = 65535.0 / rgb_max[1];
    b_inverse = 65535.0 / rgb_max[2];

    for (cur_val = c->arr; n < out_len - 1; n += 6, out += 6, cur_val += 3) {
        if (*cur_val >= rgb_max[0]) {
            r_scaled = 65535;
        } else if (*cur_val <= 0) {
            r_scaled = 0;
        } else {
            r_scaled = (uint16_t)floor(*cur_val * r_inverse);
        }
        *out = (r_scaled & 0xFF00) >> 8;
        *(out+1) = (r_scaled & 0x00FF);

        if (*(cur_val+1) >= rgb_max[1]) {
            g_scaled = 65535;
        } else if (*(cur_val+1) <= 0) {
            g_scaled = 0;
        } else {
            g_scaled = (uint16_t)floor(*(cur_val+1) * g_inverse);
        }
        *(out+2) = (g_scaled & 0xFF00) >> 8;
        *(out+3) = (g_scaled & 0x00FF);

        if (*(cur_val+2) >= rgb_max[2]) {
            b_scaled = 65535;
        } else if (*(cur_val+2) <= 0) {
            b_scaled = 0;
        } else {
            b_scaled = (uint16_t)floor(*(cur_val+2) * b_inverse);
        }
        *(out+4) = (b_scaled & 0xFF00) >> 8;
        *(out+5) = (b_scaled & 0x00FF);
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
    double *out;

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
        *out = ((double) tmp) / (double)max_val;
    }

    fclose(file);

    return c;
}
