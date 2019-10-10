#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>

#include <png.h>

#include "../linalg/linalg.h"
#include "../../color/color.h"
#include "../../color/rgb.h"

#include "canvas.h"


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
    Color rgb_tmp, srgb_tmp, srgb_max;
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

    color_default(srgb_max);

    srgb_max[0] = 1.0;
    srgb_max[1] = 1.0;
    srgb_max[2] = 1.0;

    r_inverse = 65535.0 / srgb_max[0];
    g_inverse = 65535.0 / srgb_max[1];
    b_inverse = 65535.0 / srgb_max[2];

    for (cur_val = c->arr; n < out_len - 1; n += 6, cur_val++) {
        color_copy(rgb_tmp, *cur_val);
        // clamping
        if (rgb_tmp[0] > 1.0) {
            rgb_tmp[0] = 1.0;
        } else if (rgb_tmp[0] < 0) {
            rgb_tmp[0] = 0.0;
        }
        if (rgb_tmp[1] > 1.0) {
            rgb_tmp[1] = 1.0;
        } else if (rgb_tmp[1] < 0) {
            rgb_tmp[1] = 0.0;
        }
        if (rgb_tmp[2] > 1.0) {
            rgb_tmp[2] = 1.0;
        } else if (rgb_tmp[2] < 0) {
            rgb_tmp[2] = 0.0;
        }
        rgb_to_srgb(rgb_tmp, srgb_tmp);

        if (srgb_tmp[0] > srgb_max[0]) {
            r_scaled = 65535;
        } else if (srgb_tmp[0] < 0) {
            r_scaled = 0;
        } else {
            r_scaled = (uint16_t)floor(srgb_tmp[0] * r_inverse);
        }
        *out++ = (r_scaled & 0xFF00) >> 8;
        *out++ = (r_scaled & 0x00FF);

        if (srgb_tmp[1] > srgb_max[1]) {
            g_scaled = 65535;
        } else if (srgb_tmp[1] < 0) {
            g_scaled = 0;
        } else {
            g_scaled = (uint16_t)floor(srgb_tmp[1] * g_inverse);
        }
        *out++ = (g_scaled & 0xFF00) >> 8;
        *out++ = (g_scaled & 0x00FF);

        if (srgb_tmp[2] > srgb_max[2]) {
            b_scaled = 65535;
        } else if (srgb_tmp[2] < 0) {
            b_scaled = 0;
        } else {
            b_scaled = (uint16_t)floor(srgb_tmp[2] * b_inverse);
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


/*

   mostly taken from https://gitlab.com/mmbl.pie/yar/blob/master/src/canvas_png.cpp

 */
int
write_png(Canvas c, const char *file_name)
{
  int code = 0, row;
  // NOTE: This code wraps around C code with setjmp.  I avoid C++
  // exceptions at all costs except at the end
  FILE *fp = NULL;
  uint16_t *buffer = NULL;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
	
  // Open file for writing (binary mode)
  fp = fopen(file_name, "wb");
  if (fp == NULL) {
    code = 1;
    goto cleanup;
  }

  // Initialize write structure
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    code = 2;
    goto cleanup;
  }

  // Initialize info structure
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    code = 2;
    goto cleanup;
  }

  // Setup Exception handling
  if (setjmp(png_jmpbuf(png_ptr))) {
    code = 2;
    goto cleanup;
  }

  png_init_io(png_ptr, fp);

  // Default to 48-bit RGB image
  png_set_IHDR(png_ptr, info_ptr, c->width, c->height,
               16, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  // Indicate sRGB storage
  png_set_sRGB(png_ptr, info_ptr, PNG_sRGB_INTENT_ABSOLUTE);

  png_text title_text;
  title_text.compression = PNG_TEXT_COMPRESSION_NONE;
  title_text.key = "Title";
  title_text.text = NULL;
  title_text.text_length = 0;
  png_set_text(png_ptr, info_ptr, &title_text, 1);
  png_write_info(png_ptr, info_ptr);


    buffer = (uint16_t *)malloc(c->width * c->height * 3 * sizeof(uint16_t));
    uint8_t *out = (uint8_t *)buffer;
    Color rgb_tmp, srgb_tmp, rgb_max, *cur_val;
    size_t n;
    rgb_max[0] = 1.0;
    rgb_max[1] = 1.0;
    rgb_max[2] = 1.0;
    double r_inverse = 65535.0 / rgb_max[0];
    double g_inverse = 65535.0 / rgb_max[1];
    double b_inverse = 65535.0 / rgb_max[2];
    uint16_t r_scaled, g_scaled, b_scaled;

    for (cur_val = c->arr, n = 0; n < c->width * c->height; ++n, cur_val++) {
        color_copy(rgb_tmp, *cur_val);
        // clamping
        if (rgb_tmp[0] > 1.0) {
            rgb_tmp[0] = 1.0;
        } else if (rgb_tmp[0] < 0) {
            rgb_tmp[0] = 0.0;
        }
        if (rgb_tmp[1] > 1.0) {
            rgb_tmp[1] = 1.0;
        } else if (rgb_tmp[1] < 0) {
            rgb_tmp[1] = 0.0;
        }
        if (rgb_tmp[2] > 1.0) {
            rgb_tmp[2] = 1.0;
        } else if (rgb_tmp[2] < 0) {
            rgb_tmp[2] = 0.0;
        }
        rgb_to_srgb(rgb_tmp, srgb_tmp);
        color_copy(rgb_tmp, srgb_tmp);

        if (rgb_tmp[0] > rgb_max[0]) {
            r_scaled = 65535;
        } else if (rgb_tmp[0] < 0) {
            r_scaled = 0;
        } else {
            r_scaled = (uint16_t)floor(rgb_tmp[0] * r_inverse);
        }
        *out++ = (r_scaled & 0xFF00) >> 8;
        *out++ = (r_scaled & 0x00FF);

        if (rgb_tmp[1] > rgb_max[1]) {
            g_scaled = 65535;
        } else if (rgb_tmp[1] < 0) {
            g_scaled = 0;
        } else {
            g_scaled = (uint16_t)floor(rgb_tmp[1] * g_inverse);
        }
        *out++ = (g_scaled & 0xFF00) >> 8;
        *out++ = (g_scaled & 0x00FF);

        if (rgb_tmp[2] > rgb_max[2]) {
            b_scaled = 65535;
        } else if (rgb_tmp[2] < 0) {
            b_scaled = 0;
        } else {
            b_scaled = (uint16_t)floor(rgb_tmp[2] * b_inverse);
        }
        *out++ = (b_scaled & 0xFF00) >> 8;
        *out++ = (b_scaled & 0x00FF);
    }

  for (row = 0; row < c->height; ++row) {
    png_write_row(png_ptr, (png_bytep) buffer + 3 * sizeof(uint16_t) * row * c->width);
  }

  // End write
  png_write_end(png_ptr, NULL);

 cleanup:
  if (fp != NULL) {
    fclose(fp);
  }
  if (info_ptr != NULL) {
    png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  }
  if (png_ptr != NULL) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
  }
  if (buffer != NULL) {
    free(buffer);
  }

  return code;
}

