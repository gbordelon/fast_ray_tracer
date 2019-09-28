#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

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
construct_ppm(Canvas c, bool use_clamping)
{
    int n, i, out_len;
    unsigned char *out, *buf, r_scaled, g_scaled, b_scaled;
    double *cur_val, r_inverse, g_inverse, b_inverse, r_max_val = 1.0, g_max_val = 1.0, b_max_val = 1.0;
    Ppm ppm;

    buf = (unsigned char *) malloc(header_len * sizeof(unsigned char));
    if (buf == NULL) {
        // error
        return NULL;
    }

    n = snprintf((char *)buf, header_len, "P6\n%zu %zu\n255\n", c->width, c->height);
    // error check n

    out_len = c->width * c->height * c->depth + n + 1;
    ppm = ppm_alloc(out_len);

    if (ppm == NULL) {
        free(buf);
        // error
    }
    out = ppm->arr;

    memcpy(out, buf, n * sizeof(unsigned char));
    out += n;

    free(buf);

    // iterate over all colors in c->arr and find the max
    if (use_clamping) {
        for (i = 0; i < c->width * c->height * c->depth; i += c->depth) {
            double r = c->arr[i];
            double g = c->arr[i+1];
            double b = c->arr[i+2];
            if (r > r_max_val) {
                r_max_val = r;
            }
            if (g > g_max_val) {
                g_max_val = g;
            }
            if (b > b_max_val) {
                b_max_val = b;
            }
        }
    }

    r_inverse = 255.0 / r_max_val;
    g_inverse = 255.0 / g_max_val;
    b_inverse = 255.0 / b_max_val;

    for (cur_val = c->arr; n < out_len - 1; n += 3, out += 3, cur_val += 3) {
        if (*cur_val >= r_max_val) {
            r_scaled = 255;
        } else if (*cur_val <= 0) {
            r_scaled = 0;
        } else {
            r_scaled = *cur_val * r_inverse;
        }
        *out = r_scaled;
        if (*(cur_val+1) >= g_max_val) {
            g_scaled = 255;
        } else if (*(cur_val+1) <= 0) {
            g_scaled = 0;
        } else {
            g_scaled = *(cur_val+1) * g_inverse;
        }
        *(out+1) = g_scaled;
        if (*(cur_val+2) >= b_max_val) {
            b_scaled = 255;
        } else if (*(cur_val+2) <= 0) {
            b_scaled = 0;
        } else {
            b_scaled = *(cur_val+2) * b_inverse;
        }
        *(out+2) = b_scaled;
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
