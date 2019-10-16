#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "../uthash/uthash.h"

#include "../../shapes/shapes.h"
#include "../../shapes/group.h"
#include "../../shapes/triangle.h"

#include "obj_loader.h"

struct shape_num_tuple {
    Shape shapes;
    size_t num;
};

typedef struct group_with_name {
    Shape group;
    char *name;
} group_with_name;

#define MAX_MATERIAL_NAME_LEN 256
typedef struct material_name_hash_table {
    char name[MAX_MATERIAL_NAME_LEN];
    int id;
    Material material;

    UT_hash_handle hh;
} *Named_material;

Named_material materials_ht = NULL;
    

void
set_material_flags(Material m)
{
    if (m != NULL) {
        m->reflective = m->Ks[0] > 0 || m->Ks[1] > 0 || m->Ks[2] > 0 || m->map_Ks != NULL;
        // for now link m->Tr and m->Tf
        if (m->Tr > 0 && (equal(m->Tf[0], 0) && equal(m->Tf[1], 0) && equal(m->Tf[1], 0))) {
            m->Tf[0] = m->Tr;
            m->Tf[1] = m->Tr;
            m->Tf[2] = m->Tr;
        } else if (equal(m->Tr, 0) && (m->Tf[0] > 0 || m->Tf[1] > 0 || m->Tf[2] > 0)) {
            m->Tr = (m->Tf[0] + m->Tf[1] + m->Tf[2]) / 3.0;
        }
    }
}

Pattern
parse_map(const char *line)
{
    char file_name[MAX_MATERIAL_NAME_LEN];
    size_t name_len;
    Canvas image;

    sscanf(line, "%*s %255s", file_name);
    name_len = strlen(file_name);

    // decide on uv map type, for refl maps only

    Pattern patterns = array_of_patterns(2);
    Pattern pattern_body = patterns + 1;
    if (access(file_name, F_OK ) == -1 ) {
        printf("file '%s' does not exist.", file_name);
        return NULL;
    }

    // decide on ppm or png file parser
    if (strcmp(file_name + name_len - 3, "ppm") == 0) {
        image = construct_canvas_from_ppm_file(file_name);
    } else if (strcmp(file_name + name_len - 3, "png") == 0) {
        read_png(&image, file_name);
    } else {
        printf("unrecognized file format for file %s\n", file_name);
        return NULL;
    }

    uv_texture_pattern(image, pattern_body);
    texture_map_pattern(pattern_body, TRIANGLE_UV_MAP, patterns);

    return patterns;
}

Named_material
parse_new_material(const char *line, int cur_id)
{
    Named_material s;

    s = (Named_material)malloc(sizeof(*s));
    sscanf(line, "%*s %255s", s->name);
    s->id = cur_id;
    s->material = material_alloc();

    return s;
}

void
store_material(Named_material s)
{
    if (s != NULL) {
        HASH_ADD_STR(materials_ht, name, s);
    }
}

void
parse_size_t(char *line, size_t *res)
{
    sscanf(line, "%*s %lu", res);
}

void
parse_double(char *line, double *res)
{
    sscanf(line, "%*s %lf", res);
}

void
parse_three_doubles(char *line, double *res)
{
    sscanf(line, "%*s %lf %lf %lf", res, res + 1, res + 2);
}

void
parse_mtl(FILE *mtl_file, void (*color_space_fn)(const Color, Color))
{
    char line[1024], *ptr;
    int ids = 0;
    Color tmp;
    Named_material s = NULL;
    Material current_material = NULL;
    Pattern p;

    while(fgets(line, sizeof(line), mtl_file)) {
        if (strlen(line) > 0) { // ignore empty lines
            ptr = line;
            while ((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++; // ignore leading spaces
            if (*ptr == '#' || *ptr == '\r' || *ptr == '\n' || *ptr == '\0') { // skip comments and empty lines
                continue;
            }
            if (strncmp(ptr, "newmtl", 6) == 0) {
                set_material_flags(current_material);
                store_material(s);
                s = parse_new_material(ptr, ids++);
                current_material = s->material;
            } else if (strncmp(ptr, "illum", 5) == 0) {
                parse_size_t(ptr, &(current_material->illum));
            } else if (strncmp(ptr, "d", 1) == 0) {
                parse_double(ptr, &(current_material->Tr));
                current_material->Tr = 1.0 - current_material->Tr;
            } else if (strncmp(ptr, "Tr", 2) == 0) {
                parse_double(ptr, &(current_material->Tr));
            } else if (strncmp(ptr, "Ni", 2) == 0) {
                parse_double(ptr, &(current_material->Ni));
            } else if (strncmp(ptr, "Ns", 2) == 0) {
                parse_double(ptr, &(current_material->Ns));
            } else if (strncmp(ptr, "Ka", 2) == 0) {
                parse_three_doubles(ptr, tmp);
                color_space_fn(tmp, current_material->Ka);
            } else if (strncmp(ptr, "Kd", 2) == 0) {
                parse_three_doubles(ptr, tmp);
                color_space_fn(tmp, current_material->Kd);
            } else if (strncmp(ptr, "Ks", 2) == 0) {
                parse_three_doubles(ptr, current_material->Ks);
            } else if (strncmp(ptr, "Tf", 2) == 0) {
                parse_three_doubles(ptr, current_material->Tf);
                current_material->Tf[0] = 1.0 - current_material->Tf[0];
                current_material->Tf[1] = 1.0 - current_material->Tf[1];
                current_material->Tf[2] = 1.0 - current_material->Tf[2];
            } else if (strncmp(ptr, "Ke", 2) == 0) {
                parse_three_doubles(ptr, current_material->Ke); // Ke is intensity, so no color_space transform
            } else if (strncmp(ptr, "noshadow", 8) == 0) {
                current_material->casts_shadow = false;
            } else if (strncmp(ptr, "map_Ka", 6) == 0) {
                p = parse_map(ptr);
                material_set_pattern(current_material, map_Ka, p);
            } else if (strncmp(ptr, "map_Kd", 6) == 0) {
                p = parse_map(ptr);
                material_set_pattern(current_material, map_Kd, p);
            } else if (strncmp(ptr, "map_bump", 8) == 0) {
                //p = parse_map(ptr);
                //material_set_pattern(current_material, map_bump, p);
/*
            } else if (strncmp(ptr, "map_Ks", 6) == 0) {
            } else if (strncmp(ptr, "map_Ns", 6) == 0) {
            } else if (strncmp(ptr, "map_d", 5) == 0) {
            } else if (strncmp(ptr, "bump", 4) == 0) {
            } else if (strncmp(ptr, "disp", 4) == 0) {
            } else if (strncmp(ptr, "decal", 5) == 0) {
*/
            } else {
                printf("Line \"%s\" not recognized while parsing .mtl file\n", line);
            }
        }
    }
    set_material_flags(current_material);
    store_material(s);
}

/*
f 3209/835/3210 3219/1227/3220 3220/1229/3221 3210/1219/3211 
   v0  v01 v02   v10  v11  v12  v20  v21  v22
   v0  v01 v02                  v10  v11  v12  v20  v21  v22
*/
struct shape_num_tuple
fan_triangulation(char *line, double *vertexes, double *textures, double *normals)
{
    struct shape_num_tuple rv;
    size_t v10, v11, v12, t10, t11, t12, n10, n11, n12;
    char *pch = strtok(line, " \t"), *pch2, *pch3;
    char *first_slash;
    Shape shape = array_of_shapes(32);
    size_t shape_capacity = 32;
    bool use_normals;
    bool use_textures;
    int scanf_rv;

    rv.num = 0;
    rv.shapes = shape;

    if (strchr(pch, '/') == NULL) {
        sscanf(pch, "%lu", &v10);
        use_normals = false;
        use_textures = false;
    } else {
        first_slash = strchr(pch, '/');
        if (*first_slash == '/' && *(first_slash+1) == '/') {
            sscanf(pch, "%lu//%lu", &v10, &n10);
            use_normals = true;
            use_textures = false;
        } else {
            scanf_rv = sscanf(pch, "%lu/%lu/%lu", &v10, &t10, &n10);
            if (scanf_rv == 1) {
                use_normals = false;
                use_textures = false;
            } else if (scanf_rv == 2) {
                use_normals = false;
                use_textures = true;
            } else if (scanf_rv == 3) {
                use_normals = true;
                use_textures = true;
            }
        }
    }

    pch2 = strtok(NULL, " \t");
    pch3 = strtok(NULL, " \t");
    while (pch2 != NULL && pch3 != NULL) {
        if (strchr(pch2, '/') == NULL) {
            sscanf(pch2, "%lu", &v11);
        } else {
            first_slash = strchr(pch2, '/');
            if (*first_slash == '/' && *(first_slash+1) == '/') {
                sscanf(pch2, "%lu//%lu", &v11, &n11);
            } else {
                sscanf(pch2, "%lu/%lu/%lu", &v11, &t11, &n11);
            }
        }
        if (strchr(pch3, '/') == NULL) {
            sscanf(pch3, "%lu", &v12);
        } else {
            first_slash = strchr(pch3, '/');
            if (*first_slash == '/'  && *(first_slash+1) == '/') {
                sscanf(pch3, "%lu//%lu", &v12, &n12);
            } else {
                sscanf(pch3, "%lu/%lu/%lu", &v12, &t12, &n12);
            }
        }

        // realloc if no room
        if (rv.num >= (shape_capacity - 1)) {
            rv.shapes = array_of_shapes_realloc(rv.shapes, 2 * shape_capacity);
            shape_capacity *= 2;
            shape = rv.shapes + rv.num;
        }

        if (use_normals) {
            smooth_triangle(shape,
                            vertexes + 4 * (v10-1), vertexes + 4 * (v11-1), vertexes + 4 * (v12-1),
                            normals + 4 * (n10-1), normals + 4 * (n11-1), normals + 4 * (n12-1));

        } else {
            triangle(shape,
                     vertexes + (v10-1) * 4, vertexes + (v11-1) * 4, vertexes + (v12-1) * 4);
        }

        if (use_textures) {
            memcpy(shape->fields.triangle.t1, textures + (t10-1) * 4, 3 * sizeof(double));
            memcpy(shape->fields.triangle.t2, textures + (t11-1) * 4, 3 * sizeof(double));
            memcpy(shape->fields.triangle.t3, textures + (t12-1) * 4, 3 * sizeof(double));
            shape->fields.triangle.t1[3] = 0.0;
            shape->fields.triangle.t2[3] = 0.0;
            shape->fields.triangle.t3[3] = 0.0;
            shape->fields.triangle.use_textures = true;
        }
        shape++;
        rv.num++;

        pch2 = pch3;
        pch3 = strtok(NULL, " \t");
    }

    return rv;
}

char *
parse_group(const char *line)
{
    char *rv = (char *)malloc(256 * sizeof(char));
    sscanf(line, "%*s %255s", rv);
    return rv;
}

void
parse_vertex(char *line, double *arr, size_t offset)
{
    sscanf(line, "%*s %lf %lf %lf", arr + 4 * offset, (arr + 4 * offset + 1), (arr + 4 * offset + 2));
    *(arr + 4 * offset + 3) = 1.0;
}

int
obj_parse_line(char *line, group_with_name *first_group, size_t *current_group_index, size_t *number_of_groups, double *head_vertexes, double *head_textures, double *head_normals, size_t *number_of_vertexes, size_t *number_of_textures, size_t *number_of_normals, void (*color_space_fn)(const Color, Color))
{
    static Material current_material = NULL;
    static bool started_a_group = false;
    char *new_group_name;
    struct shape_num_tuple children;
    char material_name[MAX_MATERIAL_NAME_LEN];
    char mtl_file_name[MAX_MATERIAL_NAME_LEN];
    Named_material s;
    int i;

    if (strlen(line) > 0) { // ignore empty lines
        if (strncmp(line, "v ", 2) == 0) {
            // parse vertex line
            parse_vertex(line, head_vertexes, *number_of_vertexes);
            (*number_of_vertexes)++;
            return 0;
        } else if (strncmp(line, "vt ", 3) == 0) {
            // parse texture line
            parse_vertex(line, head_textures, *number_of_textures);
            (*number_of_textures)++;
            return 1;
        } else if (strncmp(line, "vn ", 3) == 0) {
            // parse normal line
            parse_vertex(line, head_normals, *number_of_normals);
            (*number_of_normals)++;
            return 2;
        } else if (strncmp(line, "f ", 2) == 0) {
            // fan triangulate line and put triangles into the current group
            children = fan_triangulation(line+2, head_vertexes, head_textures, head_normals);
            if (children.num > 0) {
                if (current_material != NULL) {
                    for (i = 0; i < children.num; i++) {
                        shape_set_material(children.shapes + i, current_material);
                    }
                }
                started_a_group = true;
                group_add_children_stage((first_group + *current_group_index)->group,
                                   children.shapes,
                                   children.num);
            }
            free(children.shapes);//shape_free(children.shapes);
            return 3;
        } else if (strncmp(line, "g ", 2) == 0) {
            // change group context and create a new group if necessary
            new_group_name = parse_group(line);
            int i;
            for (i = 0; i < *number_of_groups; i++) {
                if (strcmp(new_group_name, (first_group + i)->name) == 0) {
                    break;
                }
            }
            if (started_a_group) {
                group_add_children_finish((first_group + *number_of_groups)->group);
                started_a_group = false;
            }
            if (i == *number_of_groups) {
                // need to allocate a new group
                (first_group + *number_of_groups)->group = array_of_shapes(1);
                group((first_group + *number_of_groups)->group, NULL, 0);
                (first_group + *number_of_groups)->name = new_group_name;
                // set current group index
                *current_group_index = *number_of_groups;
                (*number_of_groups)++;
            } else if (i < *number_of_groups) {
                // set current group index
                *current_group_index = i;
            }
            return 4;
        } else if (strncmp(line, "usemtl", 6) == 0) {
            sscanf(line, "%*s %255s", material_name);
            // lookup material by material_name in the ht
            HASH_FIND_STR(materials_ht, material_name, s);
            if (s != NULL) {
                current_material = s->material;
            } else {
                printf("Material %s not found.\n", material_name);
            }
        } else if (strncmp(line, "mtllib", 6) == 0) {
            sscanf(line, "%*s %255s", mtl_file_name);
            if (access("scenes/sibenik/sibenik.mtl", F_OK ) < 0) {
                printf("file %s not found.\n", mtl_file_name);
            } else {
                FILE *mtl_file = fopen(mtl_file_name, "r");
                if (mtl_file == NULL) {
                    printf("Error opening file %s", mtl_file_name);
                } else {
                    parse_mtl(mtl_file, color_space_fn);
                    fclose(mtl_file);
                }
            }
        }
    }
    return -1;
}

#define DEFAULT_VERTEX_NUM 32768
#define DEFAULT_GROUP_NUM 1024

void
construct_group_from_obj_file(const char *file_path, void (*color_space_fn)(const Color, Color), Shape result_group)
{
    char line[1024];
    int line_type;
    size_t vertex_count = 0, texture_count = 0, normal_count = 0;

    FILE *obj_file = fopen(file_path, "r");
    if (obj_file == NULL) {
        printf("Error opening file %s", file_path);
        return;
    }

    double *vertexes = (double *)malloc(4 * DEFAULT_VERTEX_NUM * sizeof(double));
    double *textures = (double *)malloc(4 * DEFAULT_VERTEX_NUM * sizeof(double));
    double *normals = (double *)malloc(4 * DEFAULT_VERTEX_NUM * sizeof(double));

    double *vertexes_itr = vertexes;
    double *textures_itr = textures;
    double *normals_itr = normals;

    size_t vertexes_size = DEFAULT_VERTEX_NUM;
    size_t textures_size = DEFAULT_VERTEX_NUM;
    size_t normals_size = DEFAULT_VERTEX_NUM;

    group_with_name *groups = (group_with_name *)malloc(DEFAULT_GROUP_NUM * sizeof(group_with_name));
    size_t groups_size = DEFAULT_GROUP_NUM;
    size_t groups_count = 0;
    size_t current_group_index = 0;

    groups->name = (char *) malloc((strlen("##default_group") + 1) * sizeof(char));
    strcpy(groups->name, "##default_group");
    groups->group = array_of_shapes(1);
    group(groups->group, NULL, 0);
    groups_count++;

    while(fgets(line, sizeof(line), obj_file)) {
        line_type = obj_parse_line(line,
                                groups,
                                &current_group_index, 
                                &groups_count,
                                vertexes,
                                textures,
                                normals,
                                &vertex_count,
                                &texture_count,
                                &normal_count,
                                color_space_fn);

        if (vertexes_size <= (vertex_count + 4)) {
            printf("allocating a new array\n");
            double *new_arr = (double *)malloc(4 * vertexes_size * 2 * sizeof(double));
            memcpy(new_arr, vertexes, vertex_count * 4 * sizeof(double));
            free(vertexes);
            vertexes = new_arr;
            vertexes_itr = new_arr + vertex_count * 4;
            vertexes_size *= 2;
        }
        if (textures_size <= (texture_count + 4)) {
            printf("allocating a new array\n");
            double *new_arr = (double *)malloc(4 * textures_size * 2 * sizeof(double));
            memcpy(new_arr, textures, texture_count * 4 * sizeof(double));
            free(textures);
            textures = new_arr;
            textures_itr = new_arr + texture_count * 4;
            textures_size *= 2;
        }
        if (normals_size <= (normal_count + 4)) {
            printf("allocating a new array\n");
            double *new_arr = (double *)malloc(4 * normals_size * 2 * sizeof(double));
            memcpy(new_arr, normals, normal_count * 4 * sizeof(double));
            free(normals);
            normals = new_arr;
            normals_itr = new_arr + normal_count * 4;
            normals_size *= 2;
        }
        if (groups_count >= (groups_size - 1)) {
            group_with_name *new_groups = (group_with_name *)malloc(2 * groups_size * sizeof(group_with_name));
            memcpy(new_groups, groups, groups_count * sizeof(group_with_name));
            free(groups);
            groups = new_groups;
            groups_size = 2 * groups_size;
        }
    }

    fclose(obj_file);

    // add other groups to default_group
    Shape all_groups = groups->group;
    int i;
    if (groups_count > 1) {
        all_groups = array_of_shapes(groups_count);
        for (i = 0; i < groups_count; i++) {
            memcpy(all_groups + i, (groups + i)->group, sizeof(struct shape));
        }
    }

    group(result_group, all_groups, groups_count);
    for (i = 0; i < groups_count; i++) {
        free((groups + i)->name);
    }
    free(all_groups); //shape_free

    free(groups);
    free(vertexes);
    free(textures);
    free(normals);

    Named_material s, tmp;
    HASH_ITER(hh, materials_ht, s, tmp) {
        HASH_DEL(materials_ht, s);
        free(s); // s->material is not freed because the shapes are still using them!
    }
}
