#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "../shapes/shapes.h"
#include "../light/light.h"
#include "../intersection/intersection.h"
#include "../shapes/sphere.h"
#include "../shapes/toroid.h"

#include "renderer.h"
#include "world.h"

World
world()
{
    World w = (World) malloc(sizeof(struct world));
    w->lights = NULL;
    w->lights_num = 0;
    w->shapes = NULL;
    w->shapes_num = 0;
    w->xs = intersections_empty(64);
    w->container.shapes = NULL;
    w->container.size = 0;
    w->photon_maps = NULL;
    w->global_config = NULL;

    return w;
}

/*
 * copy everything but parent pointer and intersections
 * for group and csg, allocate a shape array and recurse
 */
void
shape_copy(Shape s, Shape parent, Shape res)
{
    int i;
    Shape to, from;

    if (s == res) {
        printf("Trying to copy a shape to itself.\n");
        return;
    }

    *res = *s;
    res->parent = parent;
    res->material = NULL;
    shape_set_material(res, s->material);

    switch (s->type) {
    case SHAPE_PLANE:
    case SHAPE_SMOOTH_TRIANGLE:
    case SHAPE_TRIANGLE:
        res->xs = intersections_empty(1);
        break;
    case SHAPE_CUBE:
    case SHAPE_SPHERE:
        res->xs = intersections_empty(2);
        break;
    case SHAPE_CONE:
    case SHAPE_CYLINDER:
    case SHAPE_TOROID:
        res->xs = intersections_empty(4);
        break;
    case SHAPE_CSG:
        res->xs = intersections_empty(64);
        Shape lr = array_of_shapes(2);
        res->fields.csg.left = lr;
        res->fields.csg.right = lr+1;
        shape_copy(s->fields.csg.left, res, lr);
        shape_copy(s->fields.csg.right, res, lr+1);
        break;
    case SHAPE_GROUP:
        res->xs = intersections_empty(64);
        res->fields.group.children = array_of_shapes(s->fields.group.num_children);
        res->fields.group.children_xs = (Intersections *) malloc(s->fields.group.num_children * sizeof(Intersections));

        for (i = 0, to = res->fields.group.children, from = s->fields.group.children;
                i < s->fields.group.num_children;
                i++, to++, from++) {
            shape_copy(from, res, to);
        }
        break;
    default:
        printf("Unknown shape type %d\n", s->type);
        break;
    }
}

void
world_copy(const World w, World new_world)
{
    new_world->lights = w->lights;
    new_world->lights_num = w->lights_num;
    new_world->shapes_num = w->shapes_num;
    new_world->shapes = array_of_shapes(w->shapes_num);
    new_world->xs = intersections_empty(64);
    new_world->container.shapes = NULL;
    new_world->container.size = 0;
    new_world->photon_maps = w->photon_maps;
    new_world->global_config = w->global_config;

    int i;
    Shape from, to;
    for (i = 0, from = w->shapes, to = new_world->shapes;
            i < w->shapes_num;
            i++, from++, to++) {
        shape_copy(from, NULL, to);
    }
}

World
default_world()
{
    World w = world();
    Point p;
    point(-10, 10, -10, p);
    Color c = color(1.0, 1.0, 1.0);
    Light l = array_of_lights(1);
    point_light(p, c, l);
    w->lights = l;
    w->lights_num = 1;
    Shape shapes = (Shape) malloc(2 * sizeof(struct shape));

    Shape s1 = shapes;
    Shape s2 = shapes + 1;

    toroid(s1);
    sphere(s2);

    s1->material->Ka[0] = 0.8;
    s1->material->Ka[1] = 1.0;
    s1->material->Ka[2] = 0.6;
    s1->material->Kd[0] = 0.8;
    s1->material->Kd[1] = 1.0;
    s1->material->Kd[2] = 0.6;
    s1->material->Ks[0] = 0.0;
    s1->material->Ks[1] = 0.0;
    s1->material->Ks[2] = 0.0;
    s1->material->casts_shadow = true;
    s1->material->Tr = 0.0;
    s1->fields.toroid.r1 = 1.0;
    s1->fields.toroid.r2 = 0.5;

    Matrix tmp, tmp2;

    matrix_scale(0.5, 0.5, 0.5, tmp);
    shape_set_transform(s2, tmp);

    matrix_rotate_x(M_PI_4, tmp);
    matrix_scale(3,3,5, tmp2); 

    transform_chain(tmp, tmp2);
    shape_set_transform(s1, tmp);

    w->shapes = shapes;
    w->shapes_num = 2;

    return w;
}

Intersections
intersect_world(const World w, const Ray r, bool stop_after_first_hit)
{
    int i;

    w->xs->num = 0;
    for (i = 0; i < w->shapes_num; i++) {
        Intersections xs_1 = (w->shapes + i)->intersect(w->shapes + i, r, stop_after_first_hit);
        if (xs_1 == NULL || xs_1->num == 0) {
            continue;
        }

        // realloc
        if (xs_1->num + w->xs->num >= w->xs->array_len) {
            size_t new_array_size = (xs_1->num + w->xs->num) > (2 * w->xs->array_len)
                                  ? (xs_1->num + w->xs->num)
                                  : (2 * w->xs->array_len);
            intersections_realloc(w->xs, new_array_size);
        }

        // copy from xs_1 into xs + xs->num
        memcpy(w->xs->xs + w->xs->num, xs_1->xs, xs_1->num * sizeof(struct intersection));
        w->xs->num += xs_1->num;

        if (stop_after_first_hit) {
            break;
        }
    }

    if (w->xs->num > 1) {
        intersections_sort(w->xs);
    }

    return w->xs;
}
