#ifndef CANVAS
#define CANVAS

#include <stddef.h>
#include <stdbool.h>
#include <pthread.h>

#include "../../color/color.h"

typedef struct canvas {
    Color *arr;
    size_t width;
    size_t height;
    bool super_sample;
    void (*color_space_fn)(const Color, Color);
    pthread_mutex_t  write_lock;
} *Canvas;

typedef struct ppm_struct {
    unsigned char *arr;
    size_t len;
} *Ppm;

Canvas canvas_alloc(size_t width, size_t height, bool super_sample, void (*color_space_fn)(const Color, Color));

Ppm ppm_alloc(size_t len);

void canvas_free(Canvas c);
void ppm_free(Ppm p);

void canvas_write_pixels(Canvas c, int col, int row, Color *colors, size_t num);
void canvas_write_pixel(Canvas c, int col, int row, Color color);
void canvas_pixel_at(Canvas c, int col, int row, Color res);

int write_ppm_file(Canvas c, const bool use_scaling, const char *file_name);
int write_png(Canvas c, const char *file_name);

void construct_canvas_from_ppm_file(Canvas *c, const char * file_path, bool super_sample, void (*color_space_fn)(const Color, Color));
int read_png(Canvas *c, const char *filename, bool super_sample, void (*color_space_fn)(const Color, Color));

#endif
