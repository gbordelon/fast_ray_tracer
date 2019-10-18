#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sampler.h"

/*
 * some function from https://www.scratchapixel.com/code.php?id=34&origin=/lessons/3d-basic-rendering/global-illumination-path-tracing
 */
void
uniform_sample_hemisphere(const double r1, const double r2, Vector res) 
{ 
    double sin_theta = sqrt(1 - r1 * r1);
    double phi = 2 * M_PI * r2;
    double x = sin_theta * cos(phi);
    double z = sin_theta * sin(phi);
    Vector v;
    v[0] = x;
    v[1] = r1;
    v[2] = z;
    v[3] = 0;
    vector_normalize(v, res);
}

void
cosine_weighted_sample_hemisphere(const double r1, const double r2, Vector res)
{
    double theta = asin(sqrt(r1));
    double phi = 2 * M_PI * r2;
    Vector v;
    v[0] = sin(theta) * cos(phi);
    v[2] = sin(theta) * sin(phi);
    v[1] = cos(theta);
    v[3] = 0;
    vector_normalize(v, res);
}

void
create_coordinate_system(const Vector n, Vector nt, Vector nb) 
{ 
    Vector tmp;
    nt[3] = nb[3] = 0;
    if (fabs(n[0]) > fabs(n[1])) {
        tmp[0] = n[2];
        tmp[1] = 0;
        tmp[2] = -n[0];
        vector_scale(tmp, sqrt(n[0] * n[0] + n[2] * n[2]));
    } else {
        tmp[0] = 0;
        tmp[1] = -n[2];
        tmp[2] = n[1];
        vector_scale(tmp, sqrt(n[1] * n[1] + n[2] * n[2]));
    }
    vector_normalize(tmp, nt); // probably don't have to do this with the sqrt scaling above
    vector_cross(n, nt, nb);
} 

/*
size_t
subspace_size(const Sampler sampler, const size_t dimension)
{
    int i;
    size_t rv;
    for (i = 0, rv = 1; i < dimension; ++i) {
        rv *= sampler->steps_by_dimension[i];
    }
    return rv;
}
 * area light sample; unit square; ray with normal vector
 * point light sample
 * aperture sample; unit square

    // produce canonical representation
    for (v = 0; v < n; v++) {
        for (u = 0; u < m; u++) {
            u_jitter = (u + jitter_by(jitter)) / m
            v_jitter = (v + jitter_by(jitter)) / n
            canonical[2 * (v * m + u)] = (u + v_jitter) / m;
            canonical[2 * (v * m + u) + 1] = (v + u_jitter) / n;
        }
    }

    for (w = 0; w < o; w++) {
        for (v = 0; v < n; v++) {
            for (u = 0; u < m; u++) {
                u_jitter = (u + jitter_by(jitter)) / m
                v_jitter = (v + jitter_by(jitter)) / n
                w_jitter = (w + jitter_by(jitter)) / o
                canonical[3 * (w * n + v * m + u) + 0] = (u + (v_jitter + w_jitter)/2) / m;
                canonical[3 * (w * n + v * m + u) + 1] = (v + (u_jitter + w_jitter)/2) / n;
                canonical[3 * (w * n + v * m + u) + 2] = (w + (u_jitter + v_jitter)/2) / o;
            }
        }
    }

    // shuffle 1 D
    for (u = 0; u < m; u++) {
        k = u + jitter_by(jitter) * (m - u);
        double_swap(&canonical[u],
                    &canonical[k]);
    }

    // shuffle x and y within each x/y page
    for (w = 0; w < o; w++) {
        l = w + jitter * (o - w)
        for (v = 0; v < n; v++) {
            k = v + jitter_by(jitter) * (n - v);
            for (u = 0; u < m; u++) {
                double_swap(&canonical[3 * (w * n + v * m + u) + 0],
                            &canonical[3 * (l * n + k * m + u) + 0]);
            }
        }
        for (u = 0; u < m; u++) {
            k = u + jitter_by(jitter) * (m - u);
            for (v = 0; v < n; v++) { // add 2 to get here
                double_swap(&canonical[3 * (w * n + v * m + u) + 1],
                            &canonical[3 * (l * n + v * m + k) + 1]);
            }
        }
    }

    // shuffle x and z within each x/z page
    for (v = 0; v < n; v++) {
        l = v + jitter * (n - v)
        for (w = 0; w < o; w++) {  // add 3 to get here
            k = w + jitter_by(jitter) * (o - w);
            for (u = 0; u < m; u++) { // subtract 1 to get here
                double_swap(&canonical[3 * (w * n + v * m + u) + 0],
                            &canonical[3 * (k * n + l * m + u) + 0]);
            }
        }
        for (u = 0; u < m; u++) {
            k = u + jitter_by(jitter) * (m - u);
            for (w = 0; w < o; w++) { // add 3 to get here
                double_swap(&canonical[3 * (w * n + v * m + u) + 1],
                            &canonical[3 * (w * n + l * m + k) + 1]);
            }
        }
    }

    // shuffle y and z within each y/z page
    for (u = 0; u < m; u++) {
        l = u + jitter * (m - u)
        for (w = 0; w < o; w++) {
            k = w + jitter_by(jitter) * (o - w);
            for (v = 0; v < n; v++) {
                double_swap(&canonical[3 * (w * n + v * m + u) + 0],
                            &canonical[3 * (w * n + k * m + l) + 0]);
            }
        }
        for (v = 0; v < n; v++) {
            k = v + jitter_by(jitter) * (n - v);
            for (w = 0; w < o; w++) {
                double_swap(&canonical[3 * (w * n + v * m + u) + 1],
                            &canonical[3 * (k * n + v * m + l) + 1]);
            }
        }
    }

    // shuffle for x
    for (v = 0; v < n; v++) {
        k = v + jitter_by(jitter) * (n - v);
        for (u = 0; u < m; u++) {
            double_swap(&canonical[2 * (v * m + u) + 0],
                        &canonical[2 * (k * m + u) + 0]);
        }
    }
    // shuffle for y
    for (u = 0; u < m; u++) {
        k = u + jitter_by(jitter) * (m - u);
        for (v = 0; v < n; v++) {
            double_swap(&canonical[2 * (v * m + u) + 1],
                        &canonical[2 * (v * m + k) + 1]);
        }
    }


void
sampler_put_point(const Sampler sampler, const size_t *index, double *result)
{
    size_t index_of_first = 0;
    int i;
    for (i = sampler->dimensions - 1, index_of_first = 0; i > 0; --i) {
        index_of_first += index[i] * subspace_size(i);
    }
    index_of_first += index[0];
    memcpy(sampler->arr + sampler->dimensions * index_of_first, result, sampler->dimensions * sizeof(double));
}

void
swap_points(const Sampler sampler, const size_t *index_from, const size_t *index_to, double *tmp1, double *tmp2)
{
    sampler_get_point(sampler, index_to, tmp1);
    sampler_get_point(sampler, index_from, tmp2);
    sampler_put_point(sampler, index_to, tmp2);
    sampler_put_point(sampler, index_from, tmp1);
}

double
other_dimension_jitters_sum(const size_t total_dimensions, const size_t dimension, double * jitters)
{
    int i;
    double sum;
    for (i = sum = 0; i < total_dimensions; ++i) {
        if (i != dimension) {
            sum += jitters[i];
        }
    }
    return sum;
}

void
recursive_loop(Sampler sampler, size_t *index, double *result, double *jitters, size_t dimension)
{
    int i;
    if (dimension == 0) {
        for (i = 0; i < sampler->dimensions; ++i) {
            jitters[i] = ((double)index[i] + sampler->jitter()) / (double)sampler->steps_by_dimension[i];
        }
        for (i = 0; i < sampler->dimensions; ++i) {
            result[i] = ((double)index[i] + other_dimension_jitters_sum(sampler->dimensions, i, jitters) / (double)(sampler->dimensions-1)) / (double)sampler->steps_by_dimension[i];
        }
        sampler_put_point(sampler, index, result);
        return;
    }

    for (i = 0; i < sampler->steps_by_dimension[dimension-1]; i++) {
        index[dimension-1] = i;
        recursive_loop(sampler, index, dimension-1);
    }
}

void
correlated_shuffle(const Sampler sampler, size_t skip, size_t *index_from, size_t *index_to, const size_t dimension, const size_t depth)
{
    int i, j;
    if (dimension <= 0) {
        dimension += sampler->dimensions; // not sampler->dimensions, not depth
        if (dimension == skip) {
            --dimension;
        }
    }

    if (depth == 1) {
        if (dimension == skip) {
            --dimension;
        }
        if (dimension <= 0) {
            dimension += sampler->dimensions;
        }
        for (i = 0; i < sampler->steps_by_dimension[dimension-1]; i++) {
            index_from[dimension-1] = i;
            index_to[dimension-1] = i;
            // swap index_from and index_to
            swap_points(sampler, index_from, index_to, 
        }
        return;
    }

    for (i = 0; i < sampler->steps_by_dimension[dimension-1]; i++) {
        index_from[dimension-1] = i;
        index_to[dimension-1] = i * sampler->jitter() * (sampler->steps_by_dimension[dimension-1] - i);
        for (j = 0; j < depth - 1; ++j) {
            correlated_shuffle(sampler, skip, index_from, index_to, dimension-1-j, depth-1);
        }
    }
}

void
sampler_get_point(const Sampler sampler, const size_t *index, double *result)
{
    size_t index_of_first = 0;
    int i;
    for (i = sampler->dimensions - 1, index_of_first = 0; i > 0; --i) {
        index_of_first += index[i] * subspace_size(i);
    }
    index_of_first += index[0];
    memcpy(result, sampler->arr + sampler->dimensions * index_of_first, sampler->dimensions * sizeof(double));
}


void
sampler_reset_canonical(Sampler sampler)
{
    // canonical representation
    size_t *index_to = (size_t *)malloc(sampler->dimensions * sizeof(size_t));
    double *result = (double *)malloc(sampler->dimensions * sizeof(double));
    double *jitters = (double *)malloc(sampler->dimensions * sizeof(double));

    recursive_loop(sampler, index_to, result, jitters, sampler->dimensions);

    free(jitters);
    free(result);
    free(index_to);
}

void
sampler3D(const bool jitter, const size_t usteps, const size_t vsteps, const size_t wsteps, bool (*constraint_fn)(double *))
{
    double steps[3];
    steps[0] = usteps;
    steps[1] = vsteps;
    steps[2] = wsteps;

    sampler_init(jitter, 3, steps, constraint_fn, sampler);
}
*/

double
no_jitter()
{
    return 0.5;
}

void
double_swap(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

void
sampler_reset_canonical_2d(Sampler sampler)
{
    size_t n = sampler->steps_by_dimension[0];
    size_t m = sampler->steps_by_dimension[1];
    size_t idx;
    
    int i, j;

    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            idx = sampler->dimensions * (j * m + i);
            sampler->arr[idx    ] = (i + (j + sampler->jitter()) / (double)n) / (double)m;
            sampler->arr[idx + 1] = (j + (i + sampler->jitter()) / (double)m) / (double)n;
        }
    }
}

void
sampler_shuffle_2d(Sampler sampler)
{
    size_t m = sampler->steps_by_dimension[0];
    size_t n = sampler->steps_by_dimension[1];
    size_t swap_idx_a, swap_idx_b;
    
    int i, j, k;

    // shuffle for x
    for (j = 0; j < n; ++j) {
        k = j + sampler->jitter() * (n - j);
        for (i = 0; i < m; ++i) {
            swap_idx_a = sampler->dimensions * (j * m + i) + 0;
            swap_idx_b = sampler->dimensions * (k * m + i) + 0;
            double_swap(sampler->arr + swap_idx_a,
                        sampler->arr + swap_idx_b);
        }
    }
    // shuffle for y
    for (i = 0; i < m; ++i) {
        k = i + sampler->jitter() * (m - i);
        for (j = 0; j < n; ++j) {
            swap_idx_a = sampler->dimensions * (j * m + i) + 1;
            swap_idx_b = sampler->dimensions * (j * m + k) + 1;
            double_swap(sampler->arr + swap_idx_a,
                        sampler->arr + swap_idx_b);
        }
    }
}

void
sampler_reset_2d(Sampler sampler)
{
    sampler_reset_canonical_2d(sampler);
    sampler_shuffle_2d(sampler);
}

void
sampler_get_point_2d(Sampler sampler, const size_t *index, double *result)
{
    memcpy(result,
           sampler->arr + sampler->dimensions * (index[1] * sampler->steps_by_dimension[0] + index[0]),
           sampler->dimensions * sizeof(double));
}

void
sampler_init(const bool jitter, const size_t dimensions, const size_t *steps, bool (*constraint_fn)(const double *), Sampler sampler)
{
    if (jitter) {
        sampler->jitter = drand48;
    } else {
        sampler->jitter = no_jitter;
    }
    sampler->dimensions = dimensions;
    sampler->steps_by_dimension = (size_t *) malloc(dimensions * sizeof(size_t));
    memcpy(sampler->steps_by_dimension, steps, dimensions * sizeof(size_t));

    size_t array_len;
    int i;
    for (i = 0, array_len = 1; i < dimensions; ++i) {
        array_len *= steps[i];
    }
    sampler->arr = (double *)malloc(dimensions * array_len * sizeof(double));

    sampler->constraint_fn = constraint_fn;
    sampler->get_point = sampler_get_point_2d;
}

void
sampler_2d(const bool jitter, const size_t usteps, const size_t vsteps, bool (*constraint_fn)(const double *), Sampler sampler)
{
    size_t steps[2];
    steps[0] = usteps;
    steps[1] = vsteps;

    sampler_init(jitter, 2, steps, constraint_fn, sampler);
    sampler->reset = sampler_reset_2d;

    sampler_reset_2d(sampler);
}

void
sampler_free(Sampler sampler)
{
    if (sampler != NULL) {
        free(sampler->steps_by_dimension);
        free(sampler->arr);
    }
}

bool
sampler_default_constraint(const double *x)
{
    return true;
}


