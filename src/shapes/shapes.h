#ifndef SHAPES
#define SHAPES

#include <stdbool.h>

#include "../libs/linalg/linalg.h"
#include "../libs/canvas/canvas.h"

#include "../renderer/renderer.h"
#include "../renderer/ray.h"
#include "../intersection/intersection.h"
#include "../material/material.h"

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

struct bounding_box_fields {
    Point min; // Point
    Point max; // Point
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
    Intersections *children_xs;
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
    Point p1;
    Point p2;
    Point p3;
    Vector t1;
    Vector t2;
    Vector t3;
    Vector e1;
    Vector e2;
    bool use_textures;
    union {
        Vector normal; // triangle
        struct {
            Vector n1;
            Vector n2;
            Vector n3;
        } s_normals; // smooth triangle
    } u_normals;
};


typedef struct shape {
    Matrix transform;
    Matrix transform_inverse;
    bool transform_identity;

    Material material;
    struct shape *parent;
    Bounding_box bbox;
    Bounding_box bbox_inverse;
    bool bbox_valid;
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

    Intersections (*intersect)(struct shape *sh, Ray r, bool stop_after_first_hit);
    Intersections (*local_intersect)(struct shape *sh, Ray r, bool stop_after_first_hit);
    void (*normal_at)(struct shape *sh, Point world_point, Intersection hit, Vector res);
    void (*local_normal_at)(struct shape *sh, Point local_point, Intersection hit, Vector res);
    void (*normal_to_world)(struct shape *sh, Vector local_normal, Vector res);
    void (*world_to_object)(struct shape *sh, Point pt, Point res);

    void (*bounds)(struct shape *sh, Bounding_box *res);
    void (*parent_space_bounds)(struct shape *sh, Bounding_box *res);
    void (*divide)(struct shape *sh, size_t threshold);
    bool (*includes)(struct shape *a, struct shape *b);
} *Shape;


Shape array_of_shapes(size_t num);
Shape array_of_shapes_realloc(Shape ptr, size_t num);
void shape_free(Shape s);

void shape_normal_at(Shape sh, Point world_point, Intersection hit, Vector res);
void shape_normal_to_world(Shape sh, Vector local_normal, Vector res);
void shape_world_to_object(Shape sh, Point pt, Point res);
void shape_divide(Shape sh, size_t threshold);
bool shape_includes(Shape a, Shape b);

int shape_to_string(char *buf, size_t n, Shape sh);

Intersections shape_intersect(Shape sh, Ray r, bool stop_after_first_hit);

void shape_set_transform(Shape obj, const Matrix transform);

void shape_set_material(Shape obj, Material m);
void shape_set_material_recursive(Shape obj, Material m);

void shape_bounds(Shape sh, Bounding_box *res);
void shape_parent_space_bounds(Shape sh, Bounding_box *res);

void shape_recursive_parent_update(Shape sh, Shape parent);
void shape_recursive_invalidate_bounding_box(Shape sh);

#endif
