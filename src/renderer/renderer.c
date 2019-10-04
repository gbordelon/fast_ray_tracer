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


#include "../libs/linalg/linalg.h"
#include "../libs/sampler/sampler.h"
#include "../libs/thpool/thpool.h"
#include "../pattern/pattern.h"
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

#include "renderer.h"
#include "world.h"
#include "camera.h"
#include "ray.h"

void color_at(const World w, const Ray r, const size_t remaining, Color res, struct container *container);
void lighting(Material material, Shape shape, Light light, Point point, Vector eyev, Vector normalv, double shade_intensity, Color res);


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
ray_for_pixel(Camera cam, size_t px, size_t py, double xy_jitter[2], Ray res)
{
    double xoffset = ((double)px + xy_jitter[0]) * cam->pixel_size;
    double yoffset = ((double)py + xy_jitter[1]) * cam->pixel_size;
    double world_x = cam->half_width - xoffset;
    double world_y = cam->half_height - yoffset;
    Matrix inv;
    Point origin;
    Vector direction;
    Point p;
    Point pixel;
    Vector v;

    matrix_copy(cam->transform_inverse, inv);

    p[0] = world_x;
    p[1] = world_y;
    p[2] = -cam->canvas_distance;
    p[3] = 1.0;
    matrix_point_multiply(inv, p, pixel);

    p[0] = 0;
    p[1] = 0;
    sample_aperture(p, px, py, &cam->aperture);

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
pixel_multi_sample(Camera cam, World w, size_t x, size_t y, size_t usteps, size_t vsteps, Sampler sampler, Color res, struct container *container)
{
    double jitter[2];
    size_t u, v, index[2];
    double total_steps = (double)usteps * (double)vsteps;
    struct ray r;
    Color c;
    Color acc = color(0.0, 0.0, 0.0);

    sampler->reset(sampler);

    for (v = 0; v < vsteps; v++) {
        index[1] = v;
        for (u = 0; u < usteps; u++) {
            index[0] = u;
            sampler->get_point(sampler, index, jitter);
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
            ray_for_pixel(cam, x, y, jitter, &r);
            color_at(w, &r, 5, c, container);
            color_accumulate(acc, c);
        }
    }

    color_scale(acc, 1.0 / total_steps);
    color_copy(res, acc);
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
    size_t y_start;
    size_t y_end;
    size_t usteps;
    size_t vsteps;
    bool jitter;
    int core_id;
    Camera cam;
    Canvas image;
};

void
render_multi_helper(World w, void *args)
{
    Camera cam = ((struct render_args *)args)->cam;
    Canvas image = ((struct render_args *)args)->image;
    size_t y_start = ((struct render_args *)args)->y_start;
    size_t y_end = ((struct render_args *)args)->y_end;
    size_t usteps = ((struct render_args *)args)->usteps;
    size_t vsteps = ((struct render_args *)args)->vsteps;
    bool jitter = ((struct render_args *)args)->jitter;
    int core_id = ((struct render_args *)args)->core_id;

    //stick_this_thread_to_core(core_id);

    Color c = color(0.0,0.0,0.0);

    struct container container;
    container.shapes = NULL;
    container.size = 0;

    struct sampler sampler;
    sampler_2d(jitter, usteps, vsteps, sampler_default_constraint, &sampler);

    Color *buf = (Color *)malloc(cam->hsize * sizeof(Color));

    int i, j, k;
    for (j = y_start, k=1; j < y_end; ++j, ++k) {
        for (i = 0; i < cam->hsize; ++i) {
            color_default(c);
            pixel_multi_sample(cam, w, i, j, usteps, vsteps, &sampler, c, &container);
            color_copy(*(buf+i), c);
        }
        printf("Wrote row %lu\n", y_start);
        canvas_write_pixels(image, 0, j, buf, cam->hsize);
    }

    free(container.shapes);
    free(buf);
}

Canvas
render_multi(Camera cam, World w, size_t usteps, size_t vsteps, bool jitter, size_t num_threads)
{
    int i;
    Canvas image = canvas_alloc(cam->hsize, cam->vsize);
    World worlds = (World) malloc(num_threads * sizeof(struct world));
    struct render_args *args_array = (struct render_args *)malloc(cam->vsize * sizeof(struct render_args));
struct render_args *debug = args_array;
    for (i = 0; i < num_threads; ++i) {
        world_copy(w, worlds + i);
    }

    threadpool thpool = thpool_init(num_threads, worlds);

    for (i = 0; i < cam->vsize; ++i) {
        debug->cam = cam;
        debug->image = image;
        debug->y_start = i;
        debug->y_end = i + 1;
        debug->usteps = usteps;
        debug->vsteps = vsteps;
        debug->jitter = jitter;
        debug->core_id = i % 4; // 4 cores on this Mac
        thpool_add_work(thpool, render_multi_helper, debug);
        debug++;
    }

    thpool_wait(thpool);
    thpool_destroy(thpool);

    free(args_array);
    // TODO recursive world free for each world
    free(worlds);

    return image;
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
    struct sampler sampler;
    sampler_2d(jitter, usteps, vsteps, sampler_default_constraint, &sampler);


    k = 0;
    for (j = 0; j < cam->vsize; ++j) {
        for (i = 0; i < cam->hsize; ++i) {
            color_default(c);
            pixel_multi_sample(cam, w, i, j, usteps, vsteps, &sampler, c, &container);
            canvas_write_pixel(image, i, j, c);
        }
        k += 1;
        printf("Wrote %d rows out of %lu\n", k, cam->vsize);
    }

    free(container.shapes);

    return image;
}

void
color_at(const World w, const Ray r, const size_t remaining, Color res, struct container *container)
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


void
prepare_computations(Intersection i, Ray r, Intersections xs, Computations res, struct container *container)
{
    res->t = i->t;

    res->obj = i->object;

    ray_position(r, i->t, res->p);

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
        container->shapes = realloc(container->shapes, xs->num * sizeof(Shape *));
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
