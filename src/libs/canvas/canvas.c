#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>

#include "canvas.h"
#include "color.h"
#include "linalg.h"

#define header_len 32

Canvas
canvas_alloc(size_t width, size_t height)
{
    Canvas c = (Canvas) malloc(sizeof(struct canvas));
    // null check c
    c->arr = (Color *) malloc(width * height * sizeof(Color));
    // null check c->arr

    c->width = width;
    c->height = height;

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
canvas_write_pixels(Canvas c, int col, int row, Color *colors, size_t num)
{
    while (pthread_mutex_trylock(&(c->write_lock)) != 0);
    memcpy((c->arr + row * c->width + col), colors, num * sizeof(Color));
    pthread_mutex_unlock(&(c->write_lock));
}

void
canvas_write_pixel(Canvas c, int col, int row, Color color)
{
    while (pthread_mutex_trylock(&(c->write_lock)) != 0);
    color_copy(*(c->arr + row * c->width + col), color);
    pthread_mutex_unlock(&(c->write_lock));
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

    out_len = c->width * c->height * 3 * 2 + n + 1;
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
            if ((*(c->arr + i))[0] > rgb_max[0]) {
                rgb_max[0] = (*(c->arr + i))[0];
            }
            if ((*(c->arr + i))[1] > rgb_max[1]) {
                rgb_max[1] = (*(c->arr + i))[1];
            }
            if ((*(c->arr + i))[2] > rgb_max[2]) {
                rgb_max[2] = (*(c->arr + i))[2];
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

    total_rgb_count = width * height;
    for (i = 0, out = c->arr; i < total_rgb_count; i++, out++) {
        fscanf(file, "%u", &tmp);
        (*out)[0] = ((double) tmp) / (double)max_val;
        fscanf(file, "%u", &tmp);
        (*out)[1] = ((double) tmp) / (double)max_val;
        fscanf(file, "%u", &tmp);
        (*out)[2] = ((double) tmp) / (double)max_val;
        (*out)[3] = 0.0;
    }

    fclose(file);

    return c;
}
