#include <stdbool.h>
#include <stdlib.h>

#include "sampler.h"

double
jitter_by(const bool jitter)
{
    if (jitter) {
        return drand48();
    }
    return 0.5;
}

void
double_swap(double *a, double *b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}
