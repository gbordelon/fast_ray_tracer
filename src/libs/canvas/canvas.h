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
    pthread_mutex_t  write_lock;
} *Canvas;

typedef struct ppm_struct {
    unsigned char *arr;
    size_t len;
} *Ppm;

Canvas canvas_alloc(size_t width, size_t height);

Ppm ppm_alloc(size_t len);

void canvas_free(Canvas c);
void ppm_free(Ppm p);

void canvas_write_pixels(Canvas c, int col, int row, Color *colors, size_t num);
void canvas_write_pixel(Canvas c, int col, int row, Color color);
void canvas_pixel_at(Canvas c, int col, int row, Color res);

Ppm construct_ppm(Canvas, bool use_clamping);
int write_png(Canvas c, const char *file_name);

Canvas construct_canvas_from_ppm_file(const char * file_path);

#endif
