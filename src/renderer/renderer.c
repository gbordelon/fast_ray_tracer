#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <errno.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/sysctl.h>
#include <mach/mach_init.h>
#include <mach/thread_policy.h>
#include <sched.h>

#include "renderer.h"

#include "../libs/linalg/linalg.h"
#include "../shapes/shapes.h"
#include "../shapes/sphere.h"
#include "../shapes/plane.h"
#include "../shapes/cube.h"
#include "../shapes/cone.h"
#include "../shapes/cylinder.h"
#include "../shapes/triangle.h"
#include "../shapes/csg.h"
#include "../shapes/group.h"
#include "../shapes/toroid.h"

void
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, Color res);

double
jitter_by(bool jitter)
{
    if (jitter) {
        return drand48();
    }
    return 0.5;
}

void
area_light_point_on_light(Light l, size_t u, size_t v, double u_jitter, double v_jitter, Point retval)
{
    Vector uvec;
    vector_copy(uvec, l->u.area.uvec);

    Vector vvec;
    vector_copy(vvec, l->u.area.vvec);

    vector_scale(uvec, u_jitter);
    vector_scale(vvec, v_jitter);

    retval[0] = l->u.area.corner[0] + uvec[0] + vvec[0];
    retval[1] = l->u.area.corner[1] + uvec[1] + vvec[1];
    retval[2] = l->u.area.corner[2] + uvec[2] + vvec[2];
}

void
double_swap(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

#define AREA_LIGHT_CACHE_SIZE 65536
Points
area_light_surface_points(Light light)
{
    if (light->surface_points_cache == NULL) {
        int u, v, i;
        size_t n = light->u.area.vsteps;
        size_t m = light->u.area.usteps;
        size_t k;
        Point point = {0.0, 0.0, 0.0, 1.0};
        Points pts = (Points) malloc(AREA_LIGHT_CACHE_SIZE * sizeof(struct pts));
        Points itr = pts;
        static double *canonicalx = NULL;
        static double *canonicaly = NULL;

        if (canonicalx == NULL) {
            canonicalx = (double *) malloc(n * m * sizeof(double));
        }

        if (canonicaly == NULL) {
            canonicaly = (double *) malloc(n * m * sizeof(double));
        }

        for (i = 0, itr = pts; i < AREA_LIGHT_CACHE_SIZE; i++, itr++) {
            itr->points_num = light->num_samples;
            itr->points = (Point*) malloc(light->num_samples * sizeof(Point));

            // produce canonical representation
            for (v = 0; v < n; ++v) {
                for (u = 0; u < m; ++u) {
                    //canonicalx[(v * m + u)] = ((double)v + jitter_by(light->u.area.jitter)) / (double)n;
                    //canonicaly[(v * m + u)] = ((double)u + jitter_by(light->u.area.jitter)) / (double)m;
                    canonicalx[(v * m + u)] = ((double)u + ((double)v + jitter_by(light->u.area.jitter)) / (double)n) / (double)m;
                    canonicaly[(v * m + u)] = ((double)v + ((double)u + jitter_by(light->u.area.jitter)) / (double)m) / (double)n;

                    //canonicalx[(v * m + u)] = jitter_by(light->u.area.jitter);
                    //canonicaly[(v * m + u)] = jitter_by(light->u.area.jitter);

                    //canonical[2 * (v * m + u)] = u/(double)m + (v + jitter_by(light->u.area.jitter)) / (n*m);
                    //canonical[2 * (v * m + u) + 1] = v/(double)n + (u + jitter_by(light->u.area.jitter)) / (n*m);
                    //printf("%f %f\n", canonicalx[(v * m + u)], canonicaly[(v * m + u)]);
                }
                //printf("\n");
            }
            //printf("\n");

            // shuffle canonical for x
            for (v = 0; v < n; ++v) {
                k = v + jitter_by(light->u.area.jitter) * (double)(n - v);
                for (u = 0; u < m; ++u) {
                    double_swap(canonicalx + (v * m + u), canonicalx + (k * m + u));
                }
            }

            // shuffle canonical for y
            for (u = 0; u < m; ++u) {
                k = u + jitter_by(light->u.area.jitter) * (double)(m - u);
                for (v = 0; v < n; ++v) {
                    double_swap(canonicaly + (v * m + u), canonicaly + (v * m + k));
                }
            }

            for (v = 0; v < n; ++v) {
                for (u = 0; u < m; ++u) {
                    //printf("%f %f\n", canonicalx[(v * m + u)], canonicaly[(v * m + u)]);
                    area_light_point_on_light(light, u, v, m * canonicalx[(v * m + u)], n * canonicaly[(v * m + u)], point);
                    point_copy(*(itr->points + v * m + u), point);
                }
                //printf("\n");
            }
            //printf("\n");
        }
        light->surface_points_cache = pts;
        light->surface_points_cache_len = AREA_LIGHT_CACHE_SIZE;
        free(canonicalx);
        free(canonicaly);
    }

    int choice = rand() % light->surface_points_cache_len;
    return light->surface_points_cache + choice;
}

Points
point_light_surface_points(Light light)
{
    if (light->surface_points_cache == NULL) {
        Points pts = (Points) malloc(sizeof(struct pts));
        pts->points_num = 1;
        pts->points = (Point*) malloc(pts->points_num * sizeof(Point));
        point_copy(*pts->points, light->u.point.position);
        light->surface_points_cache = pts;
        light->surface_points_cache_len = 1;
    }
    return light->surface_points_cache;
}

double
area_light_intensity_at(Light light, World w, Point p)
{
    double total = 0.0;
    int i;
    Points pts = light->light_surface_points(light);

    for (i = 0; i < pts->points_num; i++) {
        if (!is_shadowed(w, *(pts->points + i), p)) {
            total += 1.0;
        }
    }

    return total / light->num_samples;
}

double
point_light_intensity_at(Light light, World w, Point p)
{
    if (is_shadowed(w, light->u.point.position, p)) {
        return 0.0;
    }
    return 1.0;
}

void
area_light(Point corner,
           Vector full_uvec,
           size_t usteps,
           Vector full_vvec,
           size_t vsteps,
           bool jitter,
           Color intensity,
           Light l)
{
    l->type = AREA_LIGHT;

    point_copy(l->u.area.corner, corner);

    vector_copy(l->u.area.uvec, full_uvec);
    vector_scale(l->u.area.uvec, 1.0 / (double) usteps);
    l->u.area.usteps = usteps;

    vector_copy(l->u.area.vvec, full_vvec);
    vector_scale(l->u.area.vvec, 1.0 / (double) vsteps);
    l->u.area.vsteps = vsteps;
    l->u.area.jitter = jitter;

    color_copy(l->intensity, intensity);
    l->num_samples = usteps * vsteps;


    l->light_surface_points = area_light_surface_points;
    l->intensity_at = area_light_intensity_at;

    // populate surface_points_cache
    l->surface_points_cache = NULL;
    area_light_surface_points(l);
}

Light
area_light_alloc(Point corner,
                 Vector full_uvec,
                 size_t usteps,
                 Vector full_vvec,
                 size_t vsteps,
                 bool jitter,
                 Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    area_light(corner, full_uvec, usteps, full_vvec, vsteps, jitter, intensity, l);
    return l;
}


void
point_light(Point p, Color intensity, Light l)
{
    l->type = POINT_LIGHT;
    color_copy(l->intensity, intensity);
    l->num_samples = 1;
    point_copy(l->u.point.position, p);
    l->light_surface_points = point_light_surface_points;
    l->intensity_at = point_light_intensity_at;


    // populate surface_points_cache
    l->surface_points_cache = NULL;
    point_light_surface_points(l);
}

Light
point_light_alloc(Point p, Color intensity)
{
    Light l = (Light) malloc(sizeof(struct light));
    // null check l
    point_light(p, intensity, l);
    return l;
}

Light
array_of_lights(size_t num)
{
    return (Light)malloc(num * sizeof(struct light));
}

Camera
camera(size_t hsize,
       size_t vsize,
       double field_of_view,
       double canvas_distance,
       struct aperture aperture,
       size_t sample_num,
       Matrix transform)
{
    Camera c = (Camera) malloc(sizeof(struct camera));
    // null check c

    c->hsize = hsize;
    c->vsize = vsize;
    c->field_of_view = field_of_view;
    c->canvas_distance = canvas_distance;
    c->sample_num = sample_num;
    c->aperture = aperture;


    double half_view = canvas_distance * tan(field_of_view * 0.5);
    double aspect = (double)hsize / (double)vsize;

    if (aspect >= 1.0) {
        c->half_width = half_view;
        c->half_height = half_view / aspect;
    } else {
        c->half_width = half_view * aspect;
        c->half_height = half_view;
    }

    camera_set_transform(c, transform);
    c->pixel_size = c->half_width * 2.0 / (double)hsize;

    return c;
}

void
camera_set_transform(Camera c, Matrix m)
{
    if (c != NULL) {
        matrix_copy(m, c->transform);
        matrix_inverse(m, c->transform_inverse);
    }
}

void
view_transform(Point fr, Point to, Vector up, Matrix res)
{
    Vector v;
    Vector forward;
    Vector upn;
    Vector left;
    Vector true_up;
    Matrix orientation;
    Matrix m;

    vector_from_points(to, fr, v);
    vector_normalize(v, forward);

    vector_normalize(up, upn);
    vector_cross(forward, upn, left);
    vector_cross(left, forward, true_up);

    matrix(left[0], left[1], left[2], 0,
           true_up[0], true_up[1], true_up[2], 0,
           -forward[0], -forward[1], -forward[2], 0,
           0, 0, 0, 1,
           orientation);

    matrix_translate(-fr[0], -fr[1], -fr[2], m);

    matrix_multiply(orientation, m, res);
}

World
world()
{
    World w = (World) malloc(sizeof(struct world));
    w->lights = NULL;
    w->lights_num = 0;
    w->shapes = NULL;
    w->shapes_num = 0;
    w->xs = intersections_empty(64);
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

    //memcpy(res, s, sizeof(struct shape));
    *res = *s;
    res->parent = parent;
    res->material = NULL;
    res->xs = NULL;
    res->bbox = NULL;
    res->bbox_inverse = NULL;
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
        res->fields.group.children_need_free = true;
        res->fields.group.children = array_of_shapes(s->fields.group.num_children);
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

World
world_copy(World w)
{
    World new_world = world();
    new_world->lights = w->lights;
    new_world->lights_num = w->lights_num;
    new_world->shapes_num = w->shapes_num;
    new_world->shapes = array_of_shapes(w->shapes_num);

    int i;
    Shape from, to;
    for (i = 0, from = w->shapes, to = new_world->shapes;
            i < w->shapes_num;
            i++, from++, to++) {
        shape_copy(from, NULL, to);
    }

    return new_world;
}

World
default_world()
{
    World w = world();
    Point p;
    point(-10, 10, -10, p);
    Color c = color(1.0, 1.0, 1.0);
    Light l = point_light_alloc(p, c);
    w->lights = l;
    w->lights_num = 1;
    Shape shapes = (Shape) malloc(2 * sizeof(struct shape));

    Shape s1 = shapes;
    Shape s2 = shapes + 1;

    toroid(s1);
    sphere(s2);

    s1->material->color[0] = 0.8;
    s1->material->color[1] = 1.0;
    s1->material->color[2] = 0.6;
    s1->material->diffuse = 0.7;
    s1->material->specular = 0.2;
    s1->material->casts_shadow = true;
    s1->material->transparency = 0.0;
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

bool
is_shadowed(World w, Point light_position, Point pt)
{
    Vector v;
    vector_from_points(light_position, pt, v); // from pt to light

    double distance = vector_magnitude(v); // distance between pt and light

    Vector direction;
    vector_normalize(v, direction);

    struct ray r;
    ray_array(pt, direction, &r);

    Intersections xs = intersect_world(w, &r);

    Intersection h = hit(xs, true);
    bool retval = h != NULL && h->t < distance;

    return retval;
}


void
circle_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double u, v;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
    } while (u * u + v * v > bounds->u.circle.r1);
}

void
cross_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double u, v;
    bool check;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        check = ((u > bounds->u.cross.x1) && (u <= bounds->u.cross.x2))
             || ((v > bounds->u.cross.y1) && (v <= bounds->u.cross.y2));
    } while (!check);
}

void
diamond_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double u, v;
    bool check;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        check = (u <= 0)
              ? (-u + bounds->u.diamond.b1 <= v) && (v < u + bounds->u.diamond.b2)
              : (0 <= *x)
              ? (u + bounds->u.diamond.b3 <= v) && (v < -u + bounds->u.diamond.b4)
              : false;
    } while (!check);
}

void
double_circle_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    double mag;
    double u, v;
    do {
        *x = jitter_by(true);
        *y = jitter_by(true);
        u = 2 * *x - 1;
        v = 2 * *y - 1;
        mag = u * u + v * v;
    } while (mag > bounds->u.double_circle.r1 || mag < bounds->u.double_circle.r2);
}

void
point_aperture_fn(double *x, double *y, struct aperture *bound)
{
    *x = 0.5;
    *x = 0.5;
}

void
square_aperture_fn(double *x, double *y, struct aperture *bounds)
{
    *x = jitter_by(true);
    *y = jitter_by(true);
}

void
random_point_by_function(double *x, double *y, struct aperture *bounds, void (*fn)(double *, double *, struct aperture *))
{
    fn(x, y, bounds);
    *x -= 0.5;
    *x -= 0.5;
}

void
ray_for_pixel(Camera cam, double px, double py, double x_jitter, double y_jitter, Ray res)
{
    double xoffset = (px + x_jitter) * cam->pixel_size;
    double yoffset = (py + y_jitter) * cam->pixel_size;
    double world_x = cam->half_width - xoffset;
    double world_y = cam->half_height - yoffset;
    Matrix inv;
    Point origin;
    Vector direction;
    Point p;
    Point pixel;
    Vector v;
    void (*aperture_fn)(double *, double *, struct aperture *);

    switch (cam->aperture.type) {
    case CIRCULAR_APERTURE:
        aperture_fn = circle_aperture_fn;
        break;
    case CROSS_APERTURE:
        aperture_fn = cross_aperture_fn;
        break;
    case DIAMOND_APERTURE:
        aperture_fn = diamond_aperture_fn;
        break;
    case DOUBLE_CIRCLE_APERTURE:
        aperture_fn = double_circle_aperture_fn;
        break;
    case SQUARE_APERTURE:
        aperture_fn = square_aperture_fn;
        break;
    case HEXAGONAL_APERTURE:
        // not yet impl
    case PENTAGONAL_APERTURE:
        // not yet impl
    case OCTAGONAL_APERTURE:
        // not yet impl
    case POINT_APERTURE:
    default:
        // no jittering
        // no need to do any focal_blur stuff
        aperture_fn = point_aperture_fn;
        break;
    }

    matrix_copy(cam->transform_inverse, inv);

    p[0] = world_x;
    p[1] = world_y;
    p[2] = -cam->canvas_distance;
    p[3] = 1.0;
    matrix_point_multiply(inv, p, pixel);

    p[0] = 0;
    p[1] = 0;
    random_point_by_function(p, p+1, &cam->aperture, aperture_fn);

    p[0] *= cam->aperture.size;
    p[1] *= cam->aperture.size;
    p[2] = 0;
    matrix_point_multiply(inv, p, origin);
    vector_from_points(pixel, origin, v);
    vector_normalize(v, direction);

    ray_array(origin, direction, res);
}

/*
 * ^  +-----+
 * |  |     |
 * |  |     |
 * |  +-----+
 * v u------>
 *
 * u ranges from x to x+1
 * v ranges from y to y+1
 * subdivide the pixel by units to make a grid
 * choose a point in each grid cell to sample
 * average the colors
 */
void
pixel_multi_sample(Camera cam, World w, double x, double y, size_t usteps, size_t vsteps, bool jitter, Color res, struct container *container)
{
    double x_offset, y_offset;
    size_t u, v;
    double total_steps = (double)usteps * (double)vsteps;
    double usteps_inv = 1.0 / (double)usteps;
    double vsteps_inv = 1.0 / (double)vsteps;
    struct ray r;
    Color c;
    Color acc = color(0.0, 0.0, 0.0);

    size_t n = vsteps;
    size_t m = usteps;
    size_t k;

    static double *canonical = NULL;
    if (canonical == NULL) {
        canonical = (double *) malloc(2 * vsteps * usteps * sizeof(double));
    }

    // produce canonical representation
    for (v = 0; v < n; v++) {
        for (u = 0; u < m; u++) {
            canonical[2 * (v * m + u)] = (u + (v + jitter_by(jitter)) / n) / m;
            canonical[2 * (v * m + u) + 1] = (v + (u + jitter_by(jitter)) / m) / n;
        }
    }

    // shuffle canonical for x
    for (v = 0; v < n; v++) {
        k = v + jitter_by(jitter) * (n - v);
        for (u = 0; u < m; u++) {
            double_swap(&canonical[2 * (v * m + u)], &canonical[2 * (k * m + u)]);
        }
    }

    // shuffle canonical for y
    for (u = 0; u < m; u++) {
        k = u + jitter_by(jitter) * (m - u);
        for (v = 0; v < n; v++) {
            double_swap(&canonical[2 * (v * m + u) + 1], &canonical[2 * (v * m + k) + 1]);
        }
    }

    for (v = 0; v < vsteps; v++) {
        for (u = 0; u < usteps; u++) {
            x_offset = ((double)u + canonical[2 * (v * m + u)]) * usteps_inv;
            y_offset = ((double)v + canonical[2 * (v * m + u) + 1]) * vsteps_inv;
            c[0] = 0;
            c[1] = 0;
            c[2] = 0;
            r.origin[0] = 0;
            r.origin[1] = 0;
            r.origin[2] = 0;
            r.origin[3] = 1;
            r.direction[0] = 0;
            r.direction[1] = 0;
            r.direction[2] = 0;
            r.direction[0] = 0;
            ray_for_pixel(cam, x, y, x_offset, y_offset, &r);
            color_at(w, &r, 5, c, container);
            color_accumulate(acc, c);
        }
    }

    color_scale(acc, 1.0 / total_steps);
    color_copy(res, acc);
}

Canvas
render(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter)
{
    int i,j,k;
    Color c;

    Canvas image = canvas_alloc(cam->hsize, cam->vsize);
    struct container container;
    container.shapes = NULL;
    container.size = 0;


    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            color_default(c);
            pixel_multi_sample(cam, w, (double)i, (double)j, usteps, vsteps, jitter, c, &container);
            canvas_write_pixel(image, i, j, c);
        }
        k += 1;
        printf("Wrote %d rows out of %lu\n", k, cam->vsize);
    }

    free(container.shapes);

    return image;
}

/*
 * https://yyshen.github.io/2015/01/18/binding_threads_to_cores_osx.html
 */
#define SYSCTL_CORE_COUNT   "machdep.cpu.core_count"

typedef struct cpu_set {
  uint32_t    count;
} cpu_set_t;

static inline void
CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

static inline void
CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

static inline int
CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

int
sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set)
{
  int32_t core_count = 0;
  size_t  len = sizeof(core_count);
  int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    printf("error while get core count %d\n", ret);
    return -1;
  }
  cpu_set->count = 0;
  for (int i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}

int
pthread_setaffinity_np(pthread_t thread, size_t cpu_size, cpu_set_t *cpu_set)
{
  thread_port_t mach_thread;
  int core = 0;

  for (core = 0; core < 8 * cpu_size; core++) {
    if (CPU_ISSET(core, cpu_set)) break;
  }
  printf("binding to core %d\n", core);
  thread_affinity_policy_data_t policy = { core };
  mach_thread = pthread_mach_thread_np(thread);
  thread_policy_set(mach_thread, THREAD_AFFINITY_POLICY,
                    (thread_policy_t)&policy, 1);
  return 0;
}

/*
 * https://stackoverflow.com/questions/1407786/how-to-set-cpu-affinity-of-a-particular-pthread
 */
int
stick_this_thread_to_core(int core_id) {
   int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
   if (core_id < 0 || core_id >= num_cores)
      return EINVAL;

   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(core_id, &cpuset);

   pthread_t current_thread = pthread_self();    
   return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

struct render_args {
    World w;
    Camera cam;
    Canvas image;
    size_t y_start;
    size_t y_end;
    size_t usteps;
    size_t vsteps;
    bool jitter;
    int core_id;
    int thread_id;
};

void *
render_multi_helper(void *args)
{
    Camera cam = ((struct render_args *)args)->cam;
    World w = ((struct render_args *)args)->w;
    Canvas image = ((struct render_args *)args)->image;
    size_t y_start = ((struct render_args *)args)->y_start;
    size_t y_end = ((struct render_args *)args)->y_end;
    size_t usteps = ((struct render_args *)args)->usteps;
    size_t vsteps = ((struct render_args *)args)->vsteps;
    bool jitter = ((struct render_args *)args)->jitter;
    int core_id = ((struct render_args *)args)->core_id;
    int thread_id = ((struct render_args *)args)->thread_id;

    stick_this_thread_to_core(core_id);

    Color c = color(0.0,0.0,0.0);

    struct container container;
    container.shapes = NULL;
    container.size = 0;

    Color *buf = (Color *)malloc(cam->hsize * sizeof(Color));

    int i, j, k;
    for (j = y_start, k=1; j < y_end; ++j, ++k) {
        for (i = 0; i < cam->hsize; ++i) {
            color_default(c);
            pixel_multi_sample(cam, w, (double)i, (double)j, usteps, vsteps, jitter, c, &container);
            //canvas_write_pixel(image, i, j, c);
            color_copy(*(buf+i), c);
        }
        printf("%d: Wrote %d rows out of %lu\n", thread_id, k, y_end - y_start);
        canvas_write_pixels(image, 0, j, buf, cam->hsize);
    }

    free(container.shapes);
    free(buf);
    // recursive world free

    return NULL;
}

/*
 * World is not thread safe because of Intersections buffers
 *
 * Use realloc instead of malloc to keep references in tact when a thread needs to
 * resize an Intersections array
 * 
 */
Canvas
render_multi(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter, size_t num_threads)
{
    int i;
    Canvas image = canvas_alloc(cam->hsize, cam->vsize);

    pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));

    struct render_args *args_array = (struct render_args *)malloc(num_threads * sizeof(struct render_args));

    for (i = 0; i < num_threads; i++) {
        (args_array + i)->w = world_copy(w);
        (args_array + i)->cam = cam;
        (args_array + i)->image = image;
        (args_array + i)->y_start = i * cam->vsize / num_threads;
        (args_array + i)->y_end = (i + 1) * cam->vsize / num_threads;
        (args_array + i)->usteps = usteps;
        (args_array + i)->vsteps = vsteps;
        (args_array + i)->jitter = jitter;
        (args_array + i)->core_id = i % 4;
        (args_array + i)->thread_id = i;
    }

    for (i = 0; i < num_threads; i++) {
        pthread_create(threads + i, NULL, *render_multi_helper, (void *)(args_array + i));
    }

    for (i = 0; i < num_threads; i++) {
        pthread_join(*(threads + i), NULL);
    }

    free(args_array);

    return image;
}

void
color_at(World w, Ray r, size_t remaining, Color res, struct container *container)
{
    Intersections xs = intersect_world(w, r);
    Intersection i = hit(xs, false);
    struct computations comps;
    Color c = color(0.0, 0.0, 0.0);

    if (i != NULL) {
        prepare_computations(i, r, xs, &comps, container);
        shade_hit(w, &comps, remaining, c, container);
    }
    color_copy(res, c);
}

int
sort_intersections_asc(const void *p, const void *q)
{
    double l = ((Intersection)p)->t;
    double r = ((Intersection)q)->t;
    if (l - r < 0) {
        return -1;
    } else if (l - r > 0) {
        return 1;
    }
    return 0;
}

int
sort_intersections_desc(const void *p, const void *q)
{
    double l = ((Intersection)p)->t;
    double r = ((Intersection)q)->t;
    if (l - r < 0) {
        return 1;
    } else if (l - r > 0) {
        return -1;
    }
    return 0;
}

void
intersections_reverse(Intersections xs)
{
    // sort xs by xs->xs->t descending
    qsort((void*)xs->xs, xs->num, sizeof(struct intersection), sort_intersections_desc);
}

void
intersections_sort(Intersections xs)
{
    // sort xs by xs->xs->t ascending
    qsort((void*)xs->xs, xs->num, sizeof(struct intersection), sort_intersections_asc);
}

Intersections
intersect_world(World w, Ray r)
{
    int i;

    w->xs->num = 0;
    for (i = 0; i < w->shapes_num; i++) {
        Intersections xs_1 = (w->shapes + i)->intersect(w->shapes + i, r);
        if (xs_1 == NULL || xs_1->num == 0) {
            continue;
        }

        // realloc
        if (xs_1->num + w->xs->num >= w->xs->array_len) {
            Intersections xs_2 = intersections_empty(2 * w->xs->array_len);
            memcpy(xs_2->xs, w->xs->xs, w->xs->array_len * sizeof(struct intersection));
            xs_2->num = w->xs->num;
            intersections_free(w->xs);
            w->xs = xs_2;
        }

        // copy from xs_1 into xs + xs->num
        memcpy(w->xs->xs + w->xs->num, xs_1->xs, xs_1->num * sizeof(struct intersection));
        w->xs->num += xs_1->num;
    }

    if (w->xs->num > 1) {
        intersections_sort(w->xs);
    }

    return w->xs;
}

void
position(Ray ray, double t, Point position)
{
    point_copy(position, ray->origin);
    position[0] += ray->direction[0] * t;
    position[1] += ray->direction[1] * t;
    position[2] += ray->direction[2] * t;
}

void
prepare_computations(Intersection i, Ray r, Intersections xs, Computations res, struct container *container)
{
    res->t = i->t;

    res->obj = i->object;

    position(r, i->t, res->p);

    i->object->normal_at(i->object, res->p, i, res->normalv);

    vector_copy(res->eyev, r->direction);
    vector_scale(res->eyev, -1.0);

    res->inside = false;

    vector_reflect(r->direction, res->normalv, res->reflectv);

    if (vector_dot(res->normalv, res->eyev) < 0) {
        res->inside = true;
        vector_scale(res->normalv, -1);
    }

    res->over_point[0] = res->p[0] + res->normalv[0] * EPSILON;
    res->over_point[1] = res->p[1] + res->normalv[1] * EPSILON;
    res->over_point[2] = res->p[2] + res->normalv[2] * EPSILON;
    res->over_point[3] = 1.0;

    res->under_point[0] = res->p[0] - res->normalv[0] * EPSILON;
    res->under_point[1] = res->p[1] - res->normalv[1] * EPSILON;
    res->under_point[2] = res->p[2] - res->normalv[2] * EPSILON;
    res->under_point[3] = 1.0;

    res->n1 = 1.0;
    res->n2 = 1.0;

    Computations c = res;

    int j, k;
    Intersection x;
    size_t container_len = 0;

    // check size of container array
    if (xs->num >= container->size) {
        container->shapes = realloc(container->shapes, xs->num * sizeof(Shape));
        container->size = xs->num;
    }

    for (j = 0, x = xs->xs; j < xs->num; x++, j++) {
        if (x == i) {// address compare should be okay.
            if (container_len == 0) {
                c->n1 = 1.0;
            } else {
                c->n1 = container->shapes[container_len-1]->material->refractive_index;
            }
        }

        for (k = 0; k < container_len; k++) {
            if (container->shapes[k] == i->object) {
                break;
            }
        }

        if (k < container_len) {
            // shift everything left one slot, overwriting container[index_of_object] first
            for (; k < container_len - 1; k++) {
                container->shapes[k] = container->shapes[k+1];
            }
            container_len--;
        } else {
            container->shapes[container_len] = i->object;
            container_len++;
        }

        if (x == i) {
            if (container_len == 0) {
                c->n2 = 1.0;
            } else {
                c->n2 = container->shapes[container_len-1]->material->refractive_index;
            }
            break;
        }
    }
}

void
reflected_color(World w, Computations comps, size_t remaining, Color res, struct container *container)
{
    if (remaining == 0 || equal(comps->obj->material->reflective, 0)) {
        color_default(res);
    } else {
        struct ray reflect_ray;
        Color c = color(0.0, 0.0, 0.0);
        ray_array(comps->over_point, comps->reflectv, &reflect_ray);
        color_at(w, &reflect_ray, remaining - 1, c, container);
        color_scale(c, comps->obj->material->reflective);
        color_accumulate(res, c);
    }
}

void
refracted_color(World w, Computations comps, size_t remaining, Color res, struct container *container)
{
    if (remaining == 0 || equal(comps->obj->material->transparency,0)) {
        color_default(res);
        return;
    }

    double n_ratio = comps->n1 / comps->n2;
    double cos_i = vector_dot(comps->eyev, comps->normalv);
    double sin2_t = n_ratio * n_ratio * (1.0 - cos_i * cos_i);

    if (sin2_t > 1.0) {
        color_default(res);
        return;
    }

    Color c = color(0.0, 0.0, 0.0);

    Vector t1;
    Vector t2;
    Vector direction;
    struct ray refracted_ray;

    double cos_t = sqrt(1.0 - sin2_t);
    vector_copy(t1, comps->normalv);
    vector_scale(t1, n_ratio * cos_i - cos_t);
    vector_copy(t2, comps->eyev);

    vector_scale(t2, n_ratio);
    vector(t1[0] - t2[0], t1[1] - t2[1], t1[2] - t2[2], direction);

    ray_array(comps->under_point, direction, &refracted_ray);
    color_at(w, &refracted_ray, remaining - 1, c, container);
    color_scale(c, comps->obj->material->transparency);
    color_accumulate(res, c);
}

double
schlick(Computations comps)
{
    double co = vector_dot(comps->eyev, comps->normalv);
    if (comps->n1 > comps->n2) {
        double n = comps->n1 / comps->n2;
        double sin2_t = n * n * (1.0 - co * co);
        if (sin2_t > 1.0) {
            return 1.0;
        }
        co = sqrt(1.0 - sin2_t);
    }

    double r0 = ((comps->n1 - comps->n2) / (comps->n1 + comps->n2));
    r0 = r0 * r0;

    return r0 + (1.0 - r0) * (1.0 - co) * (1.0 - co) * (1.0 - co) * (1.0 - co) * (1.0 - co);
}

void
shade_hit(World w, Computations comps, size_t remaining, Color res, struct container *container)
{
    Light itr;
    size_t i;
    Color c;
    Color surface = color(0.0, 0.0, 0.0);
    double intensity;

    for (i = 0, itr = w->lights; i < w->lights_num; i++, itr++) {
        color_default(c);
        intensity = itr->intensity_at(itr, w, comps->over_point);
        lighting(comps->obj->material,
                 comps->obj,
                 itr,
                 comps->over_point,
                 comps->eyev,
                 comps->normalv,
                 intensity,
                 c);

        color_accumulate(surface, c);
    }

    Color reflected = color(0.0, 0.0, 0.0);
    reflected_color(w, comps, remaining, reflected, container);

    Color refracted = color(0.0, 0.0, 0.0);
    refracted_color(w, comps, remaining, refracted, container);

    if (comps->obj->material->reflective > 0 && comps->obj->material->transparency > 0) {
        double reflectance = schlick(comps);
        color_scale(reflected, reflectance);
        color_scale(refracted, 1.0 - reflectance);
    }

    color_accumulate(surface, reflected);
    color_accumulate(surface, refracted);

    color_copy(res, surface);
}

void
lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, Color res)
{
    Color scolor;
    Color ambient;

    if (material->pattern != NULL) {
        material->pattern->pattern_at_shape(material->pattern, shape, point, scolor);
    } else {
        color_copy(scolor, material->color);
    }

    scolor[0] *= light->intensity[0];
    scolor[1] *= light->intensity[1];
    scolor[2] *= light->intensity[2];

    color_copy(ambient, scolor);

    color_scale(ambient, material->ambient);
    if (equal(shade_intensity, 0.0)) {
        color_accumulate(res, ambient);
        return;
    }

    Color diffuse;
    int i;
    Points pts = light->light_surface_points(light);
    Vector diff;
    Vector lightv;
    Vector reflectv;

    for (i = 0; i < pts->points_num; i++) {
        vector_from_points(*(pts->points + i), point, diff);
        vector_normalize(diff, lightv);
        double light_dot_normal = vector_dot(lightv, normalv);
        if (light_dot_normal >= 0) {
            // diffuse
            color_copy(diffuse, scolor);
            color_scale(diffuse, material->diffuse);
            color_scale(diffuse, light_dot_normal);
            color_accumulate(res, diffuse);

            // specular
            vector_scale(lightv, -1);
            vector_reflect(lightv, normalv, reflectv);
            vector_scale(lightv, -1);

            double reflect_dot_eye = vector_dot(reflectv, eyev);
            if (reflect_dot_eye > 0 && material->specular > 0) {
                double factor = pow(reflect_dot_eye, material->shininess);
                res[0] += light->intensity[0] * material->specular * factor;
                res[1] += light->intensity[1] * material->specular * factor;
                res[2] += light->intensity[2] * material->specular * factor;
            }
        }
    }

    double scaling = shade_intensity / (double)light->num_samples;
    color_scale(res, scaling);
    color_accumulate(res, ambient);
}
