#ifndef CANVAS
#define CANVAS

#include <stddef.h>
#include <stdbool.h>

typedef double Color[3];

typedef struct canvas {
    double *arr;
    size_t width;
    size_t height;
    size_t depth; // probably going to assume three
} *Canvas;

typedef struct ppm_struct {
    unsigned char *arr;
    size_t len;
} *Ppm;

Canvas canvas_alloc(size_t width, size_t height);

Ppm ppm_alloc(size_t len);

void canvas_free(Canvas c);
void ppm_free(Ppm p);

void canvas_write_pixel(Canvas c, int row, int col, Color color);

Ppm construct_ppm(Canvas, bool use_clamping);

// read ppm
Canvas construct_canvas_from_ppm_file(const char * file_path);

#endif
