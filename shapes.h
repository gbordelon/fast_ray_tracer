#ifndef SHAPES
#define SHAPES

#include <stdbool.h>
#include "linalg.h"
#include "canvas.h"

enum shape_enum {
    SHAPE_CONE,
    SHAPE_CUBE,
    SHAPE_CYLINDER,
    SHAPE_PLANE,
    SHAPE_SMOOTH_TRIANGLE,
    SHAPE_SPHERE,
    SHAPE_TRIANGLE,
    SHAPE_CSG,
    SHAPE_GROUP
};

enum csg_ops_enum {
    CSG_UNION,
    CSG_INTERSECT,
    CSG_DIFFERENCE
};

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

enum uv_map_enum {
    SPHERE_MAP,
    CYLINDER_MAP,
    PLANE_MAP,
    CUBE_MAP
};

enum uv_pattern_enum {
    UV_CHECKER,
    UV_ALIGN_CHECK,
    UV_TEXTURE
};

typedef struct pattern *Pattern;

typedef struct material {
    double color[3];
    double ambient;
    double diffuse;
    double specular;
    double shininess;
    double reflective;
    double transparency;
    double refractive_index;
    bool casts_shadow;
    Pattern pattern; // until patterns are fleshed out, lighting() should just use material->color
} *Material;

struct csg_fields {
    enum csg_ops_enum op;
    struct shape *left;
    struct shape *right;
};

struct group_fields {
    struct shape *children;
    size_t num_children;
    bool children_need_free;
    // partition children this can probably be private to the .c file
    // make subgroup this can probably be private to the .c file
};

struct cone_cylinder_fields {
    double minimum;
    double maximum;
    bool closed;
};

struct triangle_fields {
    double p1[4];
    double p2[4];
    double p3[4];
    double e1[4];
    double e2[4];
    union {
        double normal[4]; // triangle
        struct {
            double n1[4];
            double n2[4];
            double n3[4];
        } s_normals; // smooth triangle
    } u_normals;
};

typedef struct intersection *Intersection;
typedef struct intersections *Intersections;
typedef struct ray *Ray;

typedef struct shape {
    Matrix transform;
    Matrix transform_inverse;
    Material material;
    struct shape *parent;
    // a bounding box?

    enum shape_enum type;
    union {
        struct csg_fields csg;
        struct group_fields group;
        struct cone_cylinder_fields cone;
        struct cone_cylinder_fields cylinder;
        struct triangle_fields triangle;
    } fields;

    Intersections (*intersect)(struct shape *sh, Ray r);
    Intersections (*local_intersect)(struct shape *sh, Ray r);
    Vector (*normal_at)(struct shape *sh, Point world_point, Intersection hit);
    Vector (*local_normal_at)(struct shape *sh, Point local_point, Intersection hit);
    Vector (*normal_to_world)(struct shape *sh, Vector local_normal);
    Point (*world_to_object)(struct shape *sh, Point pt);

    void (*divide)(struct shape *sh, size_t threshold);
    bool (*includes)(struct shape *a, struct shape *b);
    // bounds
    // parent_space_bounds
} *Shape;

typedef struct ray {
    double origin[4];
    double direction[4];
} *Ray;

typedef struct intersection {
    double t;
    double u;
    double v;
    Shape object;
} *Intersection;

typedef struct intersections {
    Intersection xs;
    size_t array_len;
    size_t num;
} *Intersections;


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
    double a[3];
    double b[3];
};

/* concrete uv patterns */
struct uv_align_check_fields {
    double main[3];
    double ul[3];
    double ur[3];
    double bl[3];
    double br[3];
};

struct uv_checker_fields {
    double a[3];
    double b[3];
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

/* abstract UV patterns */
enum uv_map_type {
    CUBE_UV_MAP,
    CYLINDER_UV_MAP,
    PLANE_UV_MAP,
    SPHERE_UV_MAP
};

struct uv_map_pattern_fields {
    enum uv_map_type type;
    struct pattern *uv_faces;
};


typedef struct face_uv_retval {
    size_t face;
    double u;
    double v;
} *UVMapReturnType;

typedef UVMapReturnType (*uv_map_fn)(Point);

typedef struct pattern {
    Matrix transform;
    Matrix transform_inverse;
    enum pattern_type type;

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

    Color (*pattern_at_shape)(struct pattern *, Shape, Point);
    Color (*pattern_at)(struct pattern *, Point);
    Color (*uv_pattern_at)(struct pattern *, double, double);
    uv_map_fn uv_map;

    // pattern_set_transform
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

void texture_map_pattern(Pattern faces /* should be one */, enum uv_map_type type, Pattern res);


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

Pattern texture_map_pattern_alloc(Pattern faces /* should be one */, enum uv_map_type type);



void intersection(double t, Shape sh, Intersection x);
Intersection intersection_alloc(double t, Shape sh);

void intersection_with_uv(double t, double u, double v, Shape sh, Intersection x);
Intersection intersection_with_uv_alloc(double t, double u, double v, Shape sh);

Intersection hit(Intersections xs, bool filter_shadow_casters);
Intersections intersections_empty(size_t num);
void intersections_free(Intersections xs);


Ray ray_alloc(Point origin, Vector direction);
void ray_free(Ray r);
int ray_to_string(char *s, size_t n, Ray r);

// default functions for shapes
Intersections shape_intersect(Shape sh, Ray r);
Vector shape_normal_at(Shape sh, Point world_point, Intersection hit);
Vector shape_normal_to_world(Shape sh, Vector local_normal);
Point shape_world_to_object(Shape sh, Point pt);
void shape_divide(Shape sh, size_t threshold);
bool shape_includes(Shape a, Shape b);

int shape_to_string(char *buf, size_t n, Shape sh);

void shape_set_transform(Shape obj, Matrix transform);
void pattern_set_transform(Pattern pat, Matrix transform);

// Bounding Box

// Material
void material();
Material material_alloc();
void material_set_pattern(Material m, Pattern p);

#endif
