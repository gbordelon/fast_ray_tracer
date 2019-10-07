#ifndef SAMPLER
#define SAMPLER

#include <stdbool.h>
#include "../linalg/linalg.h"

/*
 * Generate a correlated multi-jitter sample set for a given unit space (plane, volume, etc)
 */
typedef struct sampler {
    size_t dimensions; // probably gonna be two or three
    size_t *steps_by_dimension;
    double *arr;

    double (*jitter)(void);
    bool (*constraint_fn)(const double *);
    void (*reset)(struct sampler *);
    void (*get_point)(struct sampler *, const size_t *, double *);
} *Sampler;

//void sampler3D(const bool jitter, const size_t usteps, const size_t vsteps, const size_t wsteps, bool (*constraint_fn)(double *));
void sampler_free(Sampler sampler);
void sampler_2d(const bool jitter, const size_t usteps, const size_t vsteps, bool (*constraint_fn)(const double *), Sampler sampler);

bool sampler_default_constraint(const double *);

/*
 * some function from https://www.scratchapixel.com/code.php?id=34&origin=/lessons/3d-basic-rendering/global-illumination-path-tracing
 */
void uniform_sample_hemisphere(const double r1, const double r2, Vector res);
void cosine_weighted_sample_hemisphere(const double r1, const double r2, Vector res);
void create_coordinate_system(const Vector n, Vector nt, Vector nb);
 


#endif
