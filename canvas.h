#ifndef CANVAS
#define CANVAS

#include <stddef.h>
#include <stdbool.h>
#include <pthread.h>

typedef double Color[4];

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

static const Color BLACK = {
    0.0,
    0.0,
    0.0,
    0.0
};

static const Color WHITE = {
    1.0,
    1.0,
    1.0,
    0.0
};

#define color_default(c) memcpy((c), BLACK, sizeof(Color));
#define color(r,g,b) { (r), (g), (b), 0.0 }

void color_accumulate(Color acc, Color other);
void color_scale(Color acc, double scalar);
void color_copy(Color to, const Color from);

int color_to_string(char *buf, size_t n, Color c);

Canvas canvas_alloc(size_t width, size_t height);

Ppm ppm_alloc(size_t len);

void canvas_free(Canvas c);
void ppm_free(Ppm p);

void canvas_write_pixels(Canvas c, int col, int row, Color *colors, size_t num);
void canvas_write_pixel(Canvas c, int col, int row, Color color);
void canvas_pixel_at(Canvas c, int col, int row, Color res);

Ppm construct_ppm(Canvas, bool use_clamping);

Canvas construct_canvas_from_ppm_file(const char * file_path);

#endif
