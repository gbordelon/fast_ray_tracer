#ifndef SHAPES
#define SHAPES

#include <stdbool.h>
#include "linalg.h"
#include "canvas.h"
#include "bounding_box.h"

enum shape_enum {
    SHAPE_CONE,
    SHAPE_CUBE,
    SHAPE_CYLINDER,
    SHAPE_PLANE,
    SHAPE_SMOOTH_TRIANGLE,
    SHAPE_SPHERE,
    SHAPE_TOROID,
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
    Pattern pattern;
    size_t ref_count;
} *Material;

struct bounding_box_fields {
    double min[4]; // Point
    double max[4]; // Point
};

struct csg_fields {
    enum csg_ops_enum op;
    struct shape *left;
    struct shape *right;
};

struct group_fields {
    struct shape *children;
    size_t num_children;
    size_t size_children_array;
    bool children_need_free;
};

struct cone_cylinder_fields {
    double minimum;
    double maximum;
    bool closed;
};

struct toroid_fields {
    double r1;
    double r2;
};

struct triangle_fields {
    double p1[4];
    double p2[4];
    double p3[4];
    double t1[4];
    double t2[4];
    double t3[4];
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
typedef struct bbox *Bounding_box; 

typedef struct shape {
    Matrix transform;
    Matrix transform_inverse;
    Material material;
    struct shape *parent;
    Bounding_box bbox;
    Bounding_box bbox_inverse;
    Intersections xs;

    enum shape_enum type;
    union {
        struct csg_fields csg;
        struct group_fields group;
        struct cone_cylinder_fields cone;
        struct cone_cylinder_fields cylinder;
        struct toroid_fields toroid;
        struct triangle_fields triangle;
    } fields;

    Intersections (*intersect)(struct shape *sh, Ray r);
    Intersections (*local_intersect)(struct shape *sh, Ray r);
    void (*normal_at)(struct shape *sh, Point world_point, Intersection hit, Vector res);
    void (*local_normal_at)(struct shape *sh, Point local_point, Intersection hit, Vector res);
    void (*normal_to_world)(struct shape *sh, Vector local_normal, Vector res);
    void (*world_to_object)(struct shape *sh, Point pt, Point res);

    Bounding_box (*bounds)(struct shape *sh);
    Bounding_box (*parent_space_bounds)(struct shape *sh);
    void (*divide)(struct shape *sh, size_t threshold);
    bool (*includes)(struct shape *a, struct shape *b);
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

/* abstract UV patterns correspond to functions in shpaes.c */
enum uv_map_type {
    CUBE_UV_MAP,
    CYLINDER_UV_MAP,
    PLANE_UV_MAP,
    SPHERE_UV_MAP,
    TOROID_UV_MAP
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

    void (*pattern_at_shape)(struct pattern *, Shape, Point, Color);
    void (*pattern_at)(struct pattern *, Shape, Point, Color);
    void (*uv_pattern_at)(struct pattern *, double, double, Color);
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

Shape array_of_shapes(size_t num);
void shape_free(Shape s);

void intersection(double t, Shape sh, Intersection x);
Intersection intersection_alloc(double t, Shape sh);

void intersection_with_uv(double t, double u, double v, Shape sh, Intersection x);
Intersection intersection_with_uv_alloc(double t, double u, double v, Shape sh);

Intersection hit(Intersections xs, bool filter_shadow_casters);
Intersections intersections_empty(size_t num);
void intersections_free(Intersections xs);


void ray_array(double origin[4], double direction[4], Ray ray);
Ray ray_alloc(Point origin, Vector direction);
void ray_free(Ray r);
int ray_to_string(char *s, size_t n, Ray r);

Intersections shape_intersect(Shape sh, Ray r);
void shape_normal_at(Shape sh, Point world_point, Intersection hit, Vector res);
void shape_normal_to_world(Shape sh, Vector local_normal, Vector res);
void shape_world_to_object(Shape sh, Point pt, Point res);
void shape_divide(Shape sh, size_t threshold);
bool shape_includes(Shape a, Shape b);

int shape_to_string(char *buf, size_t n, Shape sh);
int intersection_to_string(char *buf, size_t n, Intersection x);

void shape_set_transform(Shape obj, Matrix transform);
void pattern_set_transform(Pattern pat, Matrix transform);

void shape_set_material(Shape obj, Material m);
void shape_set_material_recursive(Shape obj, Material m);

void material();
Material material_alloc();
void material_set_pattern(Material m, Pattern p);
void material_free(Material m);

Bounding_box shape_bounds(Shape sh);
Bounding_box shape_parent_space_bounds(Shape sh);

#endif
