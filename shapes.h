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
    Matrix transform;
    
} *Pattern;

typedef struct material {
    double color[3];
    double ambient;
    double specular;
    double shininess;
    double reflective;
    double transparency;
    double refractive_index;
    bool casts_shadow;
    Pattern pattern;
} *Material;

struct csg_fields {
    enum csg_ops_enum op;
    struct shape *left;
    struct shape *right;
};

struct group_fields {
    struct shape *children;
    size_t num_children;
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
        double normal[4];
        struct {
            double n1[4];
            double n2[4]; // only for smooth triangle
            double n3[4]; // only for smooth triangle
        } s;
    } u;
};

typedef struct shape {
    Matrix transform;
    Material material;
    struct shape *parent;

    enum shape_enum type;
    union {
        struct csg_fields csg;
        struct group_fields group;
        struct cone_cylinder_fields cone;
        struct cone_cylinder_fields cylinder;
        struct triangle_fields triangle;
    } fields;

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

// Sphere
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

// Ray
// Intersection

#endif
