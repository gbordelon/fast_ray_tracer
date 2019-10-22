#ifndef PATTERN
#define PATTERN

#include "../libs/canvas/canvas.h"
#include "../libs/linalg/linalg.h"
#include "../color/color.h"
#include "../shapes/shapes.h"

enum pattern_enum {
    PATTERN_CHECKERS,
    PATTERN_GRADIENT,
    PATTERN_RADIAL_GRADIENT,
    PATTERN_RING,
    PATTERN_STRIPE,
    PATTERN_TEXTURE_MAP,
    PATTERN_CUBE_MAP,
    PATTERN_CYLINDER_MAP,
    PATTERN_BLENDED,
    PATTERN_NESTED,
    PATTERN_PERTURBED
};

enum pattern_type {
    /* concrete have pattern_at fn */
    CHECKER_PATTERN,
    GRADIENT_PATTERN,
    RADIAL_GRADIENT_PATTERN,
    RING_PATTERN,
    STRIPE_PATTERN,
    /* concrete UV instead have uv_pattern_at fn */
    UV_ALIGN_CHECKER_PATTERN,
    UV_CHECKER_PATTERN,
    UV_GRADIENT_PATTERN,
    UV_RADIAL_GRADIENT_PATTERN,
    UV_TEXTURE_PATTERN,
    /* abstract also have pattern_at_shape fn */
    BLENDED_PATTERN,
    NESTED_PATTERN,
    PERTURBED_PATTERN,
    /* abstract UV have a pattern_at fn which calls uv_pattern_at and uv_map */
    CUBE_MAP_PATTERN,
    CYLINDER_MAP_PATTERN,
    TEXTURE_MAP_PATTERN
};

/* concrete patterns */
struct concrete_pattern_fields {
    Color a;
    Color b;
};

/* concrete uv patterns */
struct uv_align_check_fields {
    Color main;
    Color ul;
    Color ur;
    Color bl;
    Color br;
};

struct uv_checker_fields {
    Color a;
    Color b;
    size_t width;
    size_t height;
};

struct uv_texture_fields {
    Canvas canvas;
};

struct pattern;

/* abstract patterns */
struct blended_pattern_fields {
    struct pattern *pattern1;
    struct pattern *pattern2;
};

struct nested_pattern_fields {
    struct pattern *pattern1;
    struct pattern *pattern2;
    struct pattern *pattern3;
};

struct perturbed_pattern_fields {
    struct pattern *pattern1;
    double frequency;
    double scale_factor;
    double persistence;
    size_t octaves;
    int seed;
};

/* abstract UV patterns correspond to functions in shpaes.c */
enum uv_map_type {
    CUBE_UV_MAP,
    CYLINDER_UV_MAP,
    PLANE_UV_MAP,
    SPHERE_UV_MAP,
    TOROID_UV_MAP,
    TRIANGLE_UV_MAP
};

struct uv_map_pattern_fields {
    enum uv_map_type type;
    struct pattern *uv_faces;
};


typedef struct face_uv_retval {
    size_t face;
    double u;
    double v;
} UVMapReturnType;

typedef void (*uv_map_fn)(Shape, Point, UVMapReturnType *);

typedef struct pattern {
    Matrix transform;
    Matrix transform_inverse;
    bool transform_identity;

    enum pattern_type type;
    size_t ref_count;

    union {
        struct concrete_pattern_fields concrete;
        struct uv_align_check_fields uv_align_check;
        struct uv_checker_fields uv_check;
        struct uv_texture_fields uv_texture;
        struct blended_pattern_fields blended;
        struct nested_pattern_fields nested;
        struct perturbed_pattern_fields perturbed;
        struct uv_map_pattern_fields uv_map;
    } fields;

    void (*pattern_at_shape)(struct pattern *, Shape, Point, ColorTriple);
    void (*pattern_at)(struct pattern *, Shape, Point, ColorTriple);
    void (*uv_pattern_at)(struct pattern *, double, double, ColorTriple);
    uv_map_fn uv_map;
} *Pattern;

void checker_pattern(Color a, Color b, Pattern res);
void gradient_pattern(Color a, Color b, Pattern res);
void radial_gradient_pattern(Color a, Color b, Pattern res);
void ring_pattern(Color a, Color b, Pattern res);
void stripe_pattern(Color a, Color b, Pattern res);

void uv_align_check_pattern(Color main, Color ul, Color ur, Color bl, Color br, Pattern res);
void uv_check_pattern(Color a, Color b, size_t width, size_t height, Pattern res);
void uv_texture_pattern(Canvas canvas, Pattern res);

void blended_pattern(Pattern p1, Pattern p2, Pattern res);
void nested_pattern(Pattern p1, Pattern p2, Pattern p3, Pattern res);
void perturbed_pattern(Pattern p1, double frequency, double scale_factor, double persistence, size_t octaves, int seed, Pattern res);

Pattern array_of_patterns(size_t num);
void texture_map_pattern(Pattern faces, enum uv_map_type type, Pattern res);


Pattern checker_pattern_alloc(Color a, Color b);
Pattern gradient_pattern_alloc(Color a, Color b);
Pattern radial_gradient_pattern_alloc(Color a, Color b);
Pattern ring_pattern_alloc(Color a, Color b);
Pattern stripe_pattern_alloc(Color a, Color b);

Pattern uv_align_check_pattern_alloc(Color main, Color ul, Color ur, Color bl, Color br);
Pattern uv_check_pattern_alloc(Color a, Color b, size_t width, size_t height);
Pattern uv_texture_pattern_alloc(Canvas canvas);

Pattern blended_pattern_alloc(Pattern p1, Pattern p2);
Pattern nested_pattern_alloc(Pattern p1, Pattern p2, Pattern p3);
Pattern perturbed_pattern_alloc(Pattern p1, double frequency, double scale_factor, double persistence, size_t octaves, int seed);

Pattern texture_map_pattern_alloc(Pattern faces, enum uv_map_type type);
void pattern_free(Pattern p);

void pattern_set_transform(Pattern pat, const Matrix transform);

#endif
