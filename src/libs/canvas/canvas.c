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
#include "../../color/srgb.h"

#include "canvas.h"


#define header_len 32

Canvas
canvas_alloc(size_t width, size_t height, bool super_sample, void (*color_space_fn)(const Color, Color))
{
    Canvas c = (Canvas) malloc(sizeof(struct canvas));
    c->arr = (Color *) malloc(width * height * sizeof(Color));

    c->width = width;
    c->height = height;
    c->color_space_fn = color_space_fn;
    c->super_sample = super_sample;

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
    if (c->super_sample) {
        Color tmp;
        color_default(tmp);
        color_default(res);
        col -= 1;
        row -= 1;
        int i, j;
        for (j = 0; j < 3; ++j, ++col) {
            if (col < 0) {
                col += c->width;
            }
            if (col == c->width) {
                col = 0;
            }
            for (i = 0; i < 3; ++i, ++row) {
                if (row < 0) {
                    row += c->height;
                }
                if (row == c->height) {
                    row = 0;
                }
                color_accumulate(tmp, *(c->arr + row * c->width + col));
            }
            row -= 3;
        }
        color_scale(tmp, 1.0 / 9.0);
        c->color_space_fn(tmp, res);
    }
    else {
        c->color_space_fn(*(c->arr + row * c->width + col), res);
    }
}

Ppm
construct_ppm(Canvas c, bool use_scaling)
{
    int n, i, out_len;
    unsigned char *out, *buf;
    uint16_t r_scaled, g_scaled, b_scaled;
    double r_inverse, g_inverse, b_inverse;
    Color rgb_tmp, srgb_tmp, srgb_max, rgb_max;
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

    color_default(rgb_max);
    for (i = 0, cur_val = c->arr; i < c->width * c->height; ++i, ++cur_val) {
        color_copy(rgb_tmp, *cur_val);
        if (rgb_tmp[0] > rgb_max[0]) {
            rgb_max[0] = rgb_tmp[0];
        }
        if (rgb_tmp[1] > rgb_max[1]) {
            rgb_max[1] = rgb_tmp[1];
        }
        if (rgb_tmp[2] > rgb_max[2]) {
            rgb_max[2] = rgb_tmp[2];
        }
    }
    print_color(rgb_max);

    color_default(srgb_max);

    for (i = 0, cur_val = c->arr; i < c->width * c->height; ++i, ++cur_val) {
        color_copy(rgb_tmp, *cur_val);
        rgb_tmp[0] /= rgb_max[0];
        rgb_tmp[1] /= rgb_max[1];
        rgb_tmp[2] /= rgb_max[2];

        rgb_to_srgb(rgb_tmp, srgb_tmp);
        if (srgb_tmp[0] > srgb_max[0]) {
            srgb_max[0] = srgb_tmp[0];
        }
        if (srgb_tmp[1] > srgb_max[1]) {
            srgb_max[1] = srgb_tmp[1];
        }
        if (srgb_tmp[2] > srgb_max[2]) {
            srgb_max[2] = srgb_tmp[2];
        }
    }

    if (rgb_max[0] < 1.0) {
        rgb_max[0] = 1.0;
    }
    if (rgb_max[1] < 1.0) {
        rgb_max[1] = 1.0;
    }
    if (rgb_max[2] < 1.0) {
        rgb_max[2] = 1.0;
    }

    print_color(srgb_max);

    r_inverse = 65535.0 / srgb_max[0];
    g_inverse = 65535.0 / srgb_max[1];
    b_inverse = 65535.0 / srgb_max[2];

    for (cur_val = c->arr; n < out_len - 1; n += 6, cur_val++) {
        color_copy(rgb_tmp, *cur_val);
        if (use_scaling) {
            rgb_tmp[0] /= rgb_max[0];
            rgb_tmp[1] /= rgb_max[1];
            rgb_tmp[2] /= rgb_max[2];
        } else {
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

int
write_ppm_file(Canvas c, const bool use_scaling, const char *file_path)
{
    FILE *pFile;
    size_t len = strlen(file_path);
    char *full_file_path = (char *)malloc((len + 5) * sizeof(char));

    strcpy(full_file_path, file_path);
    *(full_file_path + len) = '.';
    *(full_file_path + len + 1) = 'p';
    *(full_file_path + len + 2) = 'p';
    *(full_file_path + len + 3) = 'm';
    *(full_file_path + len + 4) = '\0';

    Ppm ppm = construct_ppm(c, use_scaling);

    pFile = fopen (full_file_path, "wb");
    fwrite (ppm->arr, sizeof(unsigned char), ppm->len, pFile);
    fclose (pFile);

    ppm_free(ppm);
    free(full_file_path);

    return 0;
}

void
construct_canvas_from_ppm_file(Canvas *c, const char * file_path, bool super_sample, void (*color_space_fn)(const Color, Color))
{
    char buf[32];
    size_t width, height, max_val;
    int i, total_rgb_count;
    unsigned int tmp;
    Color *out;

    FILE *file = fopen(file_path, "r");
    if (file == NULL) {
        printf("Error opening file %s", file_path);
        return;
    }
    fscanf(file, "%s", buf); // P6
    fscanf(file, "%zu", &width); // width
    fscanf(file, "%zu", &height); // height
    fscanf(file, "%zu", &max_val); // max val

    *c = canvas_alloc(width, height, super_sample, color_space_fn);
    if (*c == NULL) {
        fclose(file);
        return;
    }

    total_rgb_count = width * height;
    for (i = 0, out = (*c)->arr; i < total_rgb_count; i++, out++) {
        fscanf(file, "%u", &tmp);
        (*out)[0] = ((double) tmp) / (double)max_val;
        fscanf(file, "%u", &tmp);
        (*out)[1] = ((double) tmp) / (double)max_val;
        fscanf(file, "%u", &tmp);
        (*out)[2] = ((double) tmp) / (double)max_val;
        (*out)[3] = 0.0;
    }

    fclose(file);
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

    size_t len = strlen(file_name);
    char *full_file_path = (char *)malloc((len + 5) * sizeof(char));

    strcpy(full_file_path, file_name);
    *(full_file_path + len) = '.';
    *(full_file_path + len + 1) = 'p';
    *(full_file_path + len + 2) = 'n';
    *(full_file_path + len + 3) = 'g';
    *(full_file_path + len + 4) = '\0';

  // Open file for writing (binary mode)
  fp = fopen(full_file_path, "wb");
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
  if (full_file_path != NULL) {
    free(full_file_path);
  }

  return code;
}

int
read_png (Canvas *c, const char *filename, bool super_sample, void (*color_space_fn)(const Color, Color))
{
  int code = 0, i;
  size_t width, height;
  bool sRGB;
  png_bytep *row_ptr = NULL;
  png_structp png_ptr = NULL;
  png_byte color_type;
  png_byte bit_depth;
  png_infop info_ptr = NULL;
  uint8_t *buffer = NULL;
  // NOTE: This code wraps around C code with setjmp.  I avoid C++
  // exceptions at all costs except at the end
  FILE *fp = fopen(filename, "rb");
  if (fp == NULL) {
    code = 1;
    goto read_cleanup;
  }

  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    code = 2;
    goto read_cleanup;
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    code = 2;
    goto read_cleanup;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    code = 2;
    goto read_cleanup;
  }

  png_init_io(png_ptr, fp);

  png_read_info(png_ptr, info_ptr);

  color_type = png_get_color_type(png_ptr, info_ptr);
  bit_depth = png_get_bit_depth(png_ptr, info_ptr);
  // set image dimensions
  width = png_get_image_width(png_ptr, info_ptr);
  height = png_get_image_height(png_ptr, info_ptr);

  int intent;
  sRGB = png_get_sRGB(png_ptr, info_ptr, &intent) == PNG_INFO_sRGB;

  // we want to translate everything to 48-bit RGB color
  // http://www.libpng.org/pub/png/libpng-1.2.5-manual.html
  //if (bit_depth == 8) {
  //  png_set_filler(png_ptr, 0x00, PNG_FILLER_BEFORE);
  //}

  if (color_type == PNG_COLOR_TYPE_PALETTE) {
    png_set_palette_to_rgb(png_ptr);
  }

  if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
    png_set_expand_gray_1_2_4_to_8(png_ptr);
  }

  if (color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
    png_set_gray_to_rgb(png_ptr);
  }

  if (color_type & PNG_COLOR_MASK_ALPHA) {
    png_set_strip_alpha(png_ptr);
  }

  png_read_update_info(png_ptr, info_ptr);


  // use png_read_image interface so it can handle weird cases like
  // interlaced compression, etc
  row_ptr = (png_bytep *)malloc(height * sizeof(png_bytep));
  if (bit_depth == 8) {
    buffer = (uint8_t *)malloc(width * height * 3 * sizeof(uint8_t));
  } else if (bit_depth == 16) {
    buffer = (uint8_t *)malloc(width * height * 3 * sizeof(uint16_t));
  }
  for (i = 0; i < height; i++) {
    if (bit_depth == 8) {
      row_ptr[i] = (png_bytep) (buffer + 3 * i * width * sizeof(uint8_t));
    } else if (bit_depth == 16) {
      row_ptr[i] = (png_bytep) (buffer + 3 * i * width * sizeof(uint16_t));
    }
  }
  png_read_image(png_ptr, row_ptr);

  // convert buffer to Canvas
  *c = canvas_alloc(width, height, super_sample, color_space_fn);
  if (*c == NULL) {
    code = 3;
    goto read_cleanup;
  }

  Color *out;
  uint8_t *rgb;
  size_t total_rgb_count = width * height;
  for (i = 0, out = (*c)->arr, rgb = buffer; i < total_rgb_count; i++, out++, rgb += 3) {
    (*out)[3] = 0.0;
    if (bit_depth == 8) {
        (*out)[0] = ((double)rgb[0]) / 255.0;
        (*out)[1] = ((double)rgb[1]) / 255.0;
        (*out)[2] = ((double)rgb[2]) / 255.0;
    } else if (bit_depth == 16) {
        uint16_t tmp = (rgb[0] & 0xff) << 8;
        tmp |= rgb[1] & 0xff;
        (*out)[0] = ((double)tmp) / 65535.0;
        tmp = (rgb[2] & 0xff) << 8;
        tmp |= rgb[3] & 0xff;
        (*out)[1] = ((double)tmp) / 65535.0;
        tmp = (rgb[4] & 0xff) << 8;
        tmp |= rgb[5] & 0xff;
        (*out)[2] = ((double)tmp) / 65535.0;
        rgb += 3;
    }
  }

 read_cleanup:
  if (fp != NULL) {
    fclose(fp);
  }
  if (info_ptr != NULL) {
    png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  }
  if (png_ptr != NULL) {
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
  }
  if (row_ptr != NULL) {
    free(row_ptr);
  }
  if (buffer != NULL) {
    free(buffer);
  }

  return code;
}


