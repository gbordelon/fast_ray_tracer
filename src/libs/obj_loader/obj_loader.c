#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "obj_loader.h"

#include "../../shapes/shapes.h"
#include "../../shapes/group.h"
#include "../../shapes/triangle.h"

struct shape_num_tuple {
    Shape shapes;
    size_t num;
};

typedef struct group_with_name {
    Shape group;
    char *name;
} group_with_name;

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

    rv.num = 0;
    rv.shapes = shape;

    if (strchr(pch, '/') == NULL) {
        sscanf(pch, "%lu", &v10);
    } else {
        use_normals = true;
        first_slash = strchr(pch, '/');
        if (*first_slash == '/' && *(first_slash+1) == '/') {
            sscanf(pch, "%lu//%lu", &v10, &n10);
        } else {
            sscanf(pch, "%lu/%lu/%lu", &v10, &t10, &n10);
            use_textures = true;
        }
    }

    pch2 = strtok(NULL, " \t");
    pch3 = strtok(NULL, " \t");
    while (pch2 != NULL && pch3 != NULL) {
        if (strchr(pch2, '/') == NULL) {
            sscanf(pch2, "%lu", &v12);
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
            Shape new_arr = array_of_shapes(2 * shape_capacity);
            memcpy(new_arr, rv.shapes, rv.num * sizeof(struct shape));
            free(rv.shapes);//shape_free(rv.shapes);
            rv.shapes = new_arr;
            shape_capacity *= 2;
            shape = new_arr + rv.num;
        }

        if (use_normals) {
            smooth_triangle(shape,
                            vertexes + 4 * (v10-1), vertexes + 4 * (v11-1), vertexes + 4 * (v12-1),
                            normals + 4 * (n10-1), normals + 4 * (n11-1), normals + 4 * (n12-1));

        } else {
            triangle(shape,
                     vertexes + v10 * 4, vertexes + v11 * 4, vertexes + v12 * 4);
        }

        if (use_textures) {
            memcpy(shape->fields.triangle.t1, textures + t10 * 4, 4 * sizeof(double));
            memcpy(shape->fields.triangle.t2, textures + t11 * 4, 4 * sizeof(double));
            memcpy(shape->fields.triangle.t3, textures + t12 * 4, 4 * sizeof(double));
        }
        shape++;
        rv.num++;

        pch2 = pch3;
        pch3 = strtok(NULL, " \t");
    }

    return rv;
}

char *
parse_group(char *line)
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
obj_parse_line(char *line, group_with_name *first_group, size_t *current_group_index, size_t *number_of_groups, double *head_vertexes, double *head_textures, double *head_normals, size_t *number_of_vertexes, size_t *number_of_textures, size_t *number_of_normals)
{
    char *new_group_name;
    struct shape_num_tuple children;
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
                group_add_children((first_group + *current_group_index)->group,
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
        }
    }
    return -1;
}

#define DEFAULT_VERTEX_NUM 32768
#define DEFAULT_GROUP_NUM 1024

struct materials {
    Material arr;
    size_t size;
    size_t num;
};

void
parse_mtl(FILE *mtl_file)
{
    char line[1024];

    // I need a mechanism to map material name to the material object
    // djb2 hash function from http://www.cse.yorku.ca/~oz/hash.html
    while(fgets(line, sizeof(line), mtl_file)) {
        if (strlen(line) > 0) { // ignore empty lines
            if (strncmp(line, "newmtl", 6) == 0) { // line starts with "newmtl"
                // new material starting
                
            }
        }
    }
}

void
construct_group_from_obj_file(const char *file_path, bool use_mtl, Shape result_group)
{
    char line[1024];
    int line_type;
    size_t vertex_count = 0, texture_count = 0, normal_count = 0;

    FILE *obj_file = fopen(file_path, "r");
    if (obj_file == NULL) {
        printf("Error opening file %s", file_path);
        return;
    }

    size_t len = strlen(file_path);
    char *mtl_file_path = (char *)malloc((len + 1) * sizeof(char));

    strcpy(mtl_file_path, file_path);
    *(mtl_file_path + len - 3) = 'm';
    *(mtl_file_path + len - 2) = 't';
    *(mtl_file_path + len - 1) = 'l';
    *(mtl_file_path + len + 1) = '\0';


    FILE *mtl_file = fopen(mtl_file_path, "r");
    if (mtl_file == NULL) {
        fclose(obj_file);
        printf("Error opening file %s", mtl_file_path);
        return;
    }
    parse_mtl(mtl_file);

    fclose(mtl_file);

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
                                &normal_count);

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
    free(mtl_file_path);
}
