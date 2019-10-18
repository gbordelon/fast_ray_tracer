#include "../color/color.h"
#include "../pattern/pattern.h"

#include "material.h"

void
material(Material m)
{
    color_copy(m->Ka, WHITE);
    color_copy(m->Kd, WHITE);
    color_copy(m->Ks, WHITE);
    color_copy(m->Tf, BLACK);
    color_copy(m->Ke, BLACK);
    m->Ns = 200.0;
    m->Ni = 1.0;
    m->Tr = 0.0;
    m->reflective = false;
    m->illum = 0;
    m->casts_shadow = true;
    m->map_Ka = NULL;
    m->map_Kd = NULL;
    m->map_Ks = NULL;
    m->map_Ns = NULL;
    m->map_d = NULL;
    m->map_bump = NULL;
    m->map_disp = NULL;
    m->map_refl = NULL;

    m->ref_count = 0;
}

void
material_old(Material_old m)
{
    color_copy(m->color, WHITE);
    m->ambient = 0.1;
    m->diffuse = 0.9;
    m->specular = 0.9;
    m->shininess = 200.0;
    m->reflective = 0.0;
    m->transparency = 0.0;
    m->refractive_index = 1.0;
    m->casts_shadow = true;
    m->pattern = NULL;
    m->normal_pattern = NULL;
    m->ref_count = 0;
}

// TODO behave differently based on illum?
void
material_to_material_old(Material m, Material_old res)
{
    color_copy(res->color, m->Ka);
    res->ambient = 0.1;
    res->diffuse = 0.9;
    res->specular = 0.9;
    res->shininess = m->Ns;
    res->reflective = 0.0;
    res->transparency = m->Tr;
    res->refractive_index = m->Ni;
    res->casts_shadow = m->casts_shadow;
    res->pattern = NULL;
    res->normal_pattern = NULL;
}


Material
array_of_materials(size_t num)
{
    return (Material) malloc(num * sizeof(struct material));
}

Material_old
material_old_alloc()
{
    Material_old m = (Material_old) malloc(sizeof(struct material_old));
    material_old(m);
    return m;
}

Material
material_alloc()
{
    Material m = (Material) malloc(sizeof(struct material));
    material(m);
    m->ref_count++;
    return m;
}

void
material_free(Material m)
{
    if (m != NULL) {
        m->ref_count--;
        if (m->ref_count == 0) {
            if (m->map_Ka != NULL) {
                pattern_free(m->map_Ka);
            }
            if (m->map_Kd != NULL) {
                pattern_free(m->map_Kd);
            }
            if (m->map_Ks != NULL) {
                pattern_free(m->map_Ks);
            }
            if (m->map_Ns != NULL) {
                pattern_free(m->map_Ns);
            }
            if (m->map_d != NULL) {
                pattern_free(m->map_d);
            }
            if (m->map_bump != NULL) {
                pattern_free(m->map_bump);
            }
            if (m->map_disp != NULL) {
                pattern_free(m->map_disp);
            }
            if (m->map_refl != NULL) {
                pattern_free(m->map_refl);
            }

            free(m);
        }
    }
}
