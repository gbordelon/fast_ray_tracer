#ifndef SHAPES
#define SHAPES

#include <stdbool.h>
#include "linalg.h"

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

typedef struct pattern {
    Matrix transform; // struct to include transform and its inverse
    // pattern_at
    // pattern_at_shape
    // what else? probably look like Shape
} *Pattern;

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
    // add child
    // add children
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
    Matrix transform; // maybe a struct to hold a transform and its inverse
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
    Point (*world_to_object)(Point pt);

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
    size_t len;
} *Intersections;

// intersect
// local_intersect for each shape

Intersection hit(Intersections xs);
Intersections intersections_empty(size_t num);
void intersections_free(Intersections xs);


Ray ray_alloc(Point origin, Vector direction);

// Sphere
Shape sphere();

// Plane
// Cube
// Cone
// Cylinder
// Triangle
// Smooth Triangle
// CSG
// Group

// Bounding Box
// Material
// Pattern
// SphereUVMap
// CylinderUVMap
// PlaneUVMap
// CubeUVMap
// TextureMapPattern
// CubeMapPattern
// CylinderMapPattern
// BlendedPattern
// NestedPattern
// PerturbedPattern
// UVCheckerPattern
// UVAlignCheckPattern
// UVTexturePattern
// Checkers
// Gradient
// RadialGradient
// Ring
// Stripe

#endif
