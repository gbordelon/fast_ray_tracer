#include <math.h>

#include "color.h"
#include "lab.h"
#include "xyz.h"

void
lab_to_xyz(const Color lab, Color xyz)
{
    double p = (lab[0] + 16.0) / 116.0;
    xyz[0] = tristimulus[0] * pow(p + lab[1] / 500.0, 3.0);
    xyz[1] = tristimulus[1] * pow(p, 3.0);
    xyz[2] = tristimulus[2] * pow(p - lab[2] / 200.0, 3.0);
}

void
lab_to_rgb(const Color lab, Color rgb)
{
    Color xyz;
    lab_to_xyz(lab, xyz);
    xyz_to_rgb(xyz, rgb);
}

int
lab_compare_l(const Color l, const Color r)
{
    if (l[0] - r[0] < 0) {
        return -1;
    } else if (l[0] - r[0] > 0) {
        return 1;
    }
    return 0;
}

int
lab_compare_a(const Color l, const Color r)
{
    if (l[1] - r[1] < 0) {
        return -1;
    } else if (l[1] - r[1] > 0) {
        return 1;
    }
    return 0;
}

int
lab_compare_b(const Color l, const Color r)
{
    if (l[2] - r[2] < 0) {
        return -1;
    } else if (l[2] - r[2] > 0) {
        return 1;
    }
    return 0;
}
