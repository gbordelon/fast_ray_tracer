#ifndef CAMERA
#define CAMERA

#include <stdbool.h>

#include "../libs/linalg/linalg.h"
#include "../libs/sampler/sampler.h"

enum aperture_type {
    CIRCULAR_APERTURE,
    CROSS_APERTURE,
    DIAMOND_APERTURE,
    DOUGHNUT_APERTURE,
    HEXAGONAL_APERTURE,
    PENTAGONAL_APERTURE,
    POINT_APERTURE,
    SQUARE_APERTURE,
    OCTAGONAL_APERTURE,
};

struct circle_aperture_args {
    double r1;
};

struct cross_aperture_args {
    double x1;
    double x2;
    double y1;
    double y2;
};

struct diamond_aperture_args {
    double b1;
    double b2;
    double b3;
    double b4;
};

struct doughnut_aperture_args {
    double r1;
    double r2;
};

typedef struct aperture {
    enum aperture_type type;
    double size;
    bool jitter;
    struct sampler sampler;

    union {
        struct circle_aperture_args circle;
        struct cross_aperture_args cross;
        struct diamond_aperture_args diamond;
        struct doughnut_aperture_args doughnut;
    } u;

    void (*aperture_fn)(double *, double *, struct aperture *);
} *Aperture;

typedef struct camera {
    size_t hsize;
    size_t vsize;
    size_t usteps;
    size_t vsteps;
    double field_of_view;
    double canvas_distance;
    struct aperture aperture;
    double half_width;
    double half_height;
    double pixel_size;
    Matrix transform;
    Matrix transform_inverse;
} *Camera;

Camera camera(size_t hsize, size_t vsize, double field_of_view, double canvas_distance, size_t usteps, size_t vsteps, Aperture aperture, Matrix transform);
void view_transform(Point fr, Point to, Vector up, Matrix res);
void camera_set_transform(Camera c, Matrix m);
void aperture(enum aperture_type type, double size, size_t usteps, size_t vsteps, bool jitter, Aperture res);

void sample_aperture(double xy[2], size_t u, size_t v, const Aperture aperture);

void circle_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct circle_aperture_args *args, Aperture res);
void cross_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct cross_aperture_args *args, Aperture res);
void diamond_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct diamond_aperture_args *args, Aperture res);
void doughnut_aperture(double size, size_t usteps, size_t vsteps, bool jitter, struct doughnut_aperture_args *args, Aperture res);
void square_aperture(double size, size_t usteps, size_t vsteps, bool jitter, Aperture res);
void point_aperture(Aperture res);

#endif
