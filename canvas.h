#ifndef CANVAS
#define CANVAS

#include <stddef.h>
#include <stdbool.h>

typedef struct color {
    double arr[3];
} *Color;

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

Color color_default();
Color color(double r, double g, double b);
void color_accumulate(Color acc, Color other);
void color_scale(Color acc, double scalar);
int color_to_string(char *buf, size_t n, Color c);

Canvas canvas_alloc(size_t width, size_t height);

Ppm ppm_alloc(size_t len);

void canvas_free(Canvas c);
void ppm_free(Ppm p);
void color_free(Color c);

void canvas_write_pixel(Canvas c, int col, int row, Color color);
void canvas_pixel_at(Canvas c, int col, int row, Color res);
Color canvas_pixel_at_alloc(Canvas c, int col, int row);

Ppm construct_ppm(Canvas, bool use_clamping);

// read ppm
Canvas construct_canvas_from_ppm_file(const char * file_path);

#endif
