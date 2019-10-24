#ifndef SAMPLER
#define SAMPLER

#include <stdbool.h>
#include "../linalg/linalg.h"

/*
 * Generate a correlated multi-jitter sample set for a given unit space (plane, volume, etc)
 */
typedef struct sampler {
    size_t dimensions; // probably gonna be two or three
    bool needs_hemi_coords;
    Vector nt, nb;
    size_t *steps_by_dimension;
    double *arr;

    double (*jitter)(void);
    bool (*constraint_fn)(const double *);
    void (*reset)(struct sampler *);
    void (*get_point)(struct sampler *, const size_t *, double *);
    void (*get_vector_hemisphere)(struct sampler *, Vector, bool, size_t *, double *, Vector);
    void (*reset_hemisphere_coords)(struct sampler *);
} *Sampler;

//void sampler3D(const bool jitter, const size_t usteps, const size_t vsteps, const size_t wsteps, bool (*constraint_fn)(double *));
void sampler_free(Sampler sampler);
void sampler_2d(const bool jitter, const size_t usteps, const size_t vsteps, bool (*constraint_fn)(const double *), Sampler sampler);

bool sampler_default_constraint(const double *);

#endif
