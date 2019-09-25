#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "obj_loader.h"
#include "shapes.h"
#include "group.h"
#include "triangle.h"

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
fan_triangulation(char *line, double *vertexes, double *textures, double *normals, bool contains_textures, bool contains_normals)
{
    struct shape_num_tuple rv;
    size_t v10, v11, v12, t10, t11, t12, n10, n11, n12;
    char *pch = strtok(NULL, " \t"), *pch2, *pch3;
    char *first_slash;
    Shape shape = array_of_shapes(32);
    size_t shape_capacity = 32;

    rv.num = 0;
    rv.shapes = shape;

    if (strchr(pch, '/') == NULL) {
        sscanf(pch, "%lu", &v10);
    } else {
        first_slash = strchr(pch, '/');
        if (*first_slash == '/' && *(first_slash+1) == '/') {
            sscanf(pch, "%lu//%lu", &v10, &n10);
        } else {
            sscanf(pch, "%lu/%lu/%lu", &v10, &t10, &n10);
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
        if (rv.num == shape_capacity) {
            Shape new_arr = array_of_shapes(2 * shape_capacity);
            memcpy(new_arr, rv.shapes, rv.num * sizeof(struct shape));
            shape_free(rv.shapes);
            rv.shapes = new_arr;
            shape_capacity *= 2;
            shape = new_arr + rv.num;
        }

        if (contains_normals) {
            smooth_triangle(shape,
                            vertexes + v10 * 4, vertexes + v11 * 4, vertexes + v12 * 4,
                            normals + v10 * 4, vertexes + v11 * 4, vertexes + v12 * 4);

        } else {
            triangle(shape,
                     vertexes + v10 * 4, vertexes + v11 * 4, vertexes + v12 * 4);
        }

        if (contains_textures) {
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
parse_vertex(char *line, double **arr)
{
    sscanf(line, "%*s %lf %lf %lf", *arr, *(arr + 1), *(arr + 2));
    **(arr + 3) = 1.0;
    (*arr) += 4;
}

int
parse_line(char *line, group_with_name *first_group, size_t *current_group_index, size_t *number_of_groups, double **vertexes, double **textures, double **normals, double *head_vertexes, double *head_textures, double *head_normals, size_t *number_of_vertexes, size_t *number_of_textures, size_t *number_of_normals)
{
    char *pch, *new_group_name;
    struct shape_num_tuple children;
    pch = strtok(line, " \t");
    if (pch != NULL) { // ignore empty lines
        switch (*pch) {
        case 'v':
            if (*(pch+1) == 't') {
                // parse texture line
                parse_vertex(line, textures);
                (*number_of_textures)++;
                return 1;
            } else if (*(pch+1) == 'n') {
                // parse normal line
                parse_vertex(line, normals);
                (*number_of_normals)++;
                return 2;
            }
            // parse vertex line
            parse_vertex(line, vertexes);
            (*number_of_vertexes)++;
            return 0;
            break;
        case 'f':
            // fan triangulate line and put triangles into the current group
            children = fan_triangulation(line+2, head_vertexes, head_textures, head_normals, *number_of_textures > 0, *number_of_normals > 0);
            group_add_children((first_group + *current_group_index)->group,
                               children.shapes,
                               children.num);
            return 3;
            break;
        case 'g':
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
                (first_group + *number_of_groups)->name = new_group_name;
                // set current group index
                *current_group_index = *number_of_groups;
                (*number_of_groups)++;
            } else if (i < *number_of_groups) {
                // set current group index
                *current_group_index =i;
            }
            return 4;
            break;
        default:
            // ignore lines starting with #
            return -2;
            break;
        }
    }
    return -1;
}

#define DEFAULT_VERTEX_NUM 1024
#define DEFAULT_GROUP_NUM 1024

Shape
construct_group_from_obj_file(const char *file_path)
{
    char line[1024];
    int line_type;
    size_t vertex_count = 0, texture_count = 0, normal_count = 0;

    FILE *file = fopen(file_path, "r");
    if (file == NULL) {
        printf("Error opening file %s", file_path);
        return NULL;
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

    groups->group = array_of_shapes(1);
    groups->name = (char *) malloc((strlen("##default_group") + 1) * sizeof(char));
    strcpy(groups->name, "##default_group");
    group(groups->group, NULL, 0);
    groups_count++;

    while(fgets(line, sizeof(line), file)) {
        line_type = parse_line(line,
                                groups,
                                &current_group_index, 
                                &groups_count,
                                &vertexes_itr,
                                &textures_itr,
                                &normals_itr,
                                vertexes,
                                textures,
                                normals,
                                &vertex_count,
                                &texture_count,
                                &normal_count);
        switch (line_type) {
        case 0:
            break;
        case 1:
            break;;
        case 2:
            break;
        case 3:
            break;
        case 4:
            break;
        default:
            break;
        }
        if (vertexes_size <= (vertex_count + 4)) {
            double *new_arr = (double *)malloc(4 * vertexes_size * 2 * sizeof(double));
            memcpy(new_arr, vertexes, vertex_count * 4 * sizeof(double));
            free(vertexes);
            vertexes = new_arr;
            vertexes_itr = new_arr + vertex_count * 4;
            vertexes_size *= 2;
        }
        if (textures_size <= (texture_count + 4)) {
            double *new_arr = (double *)malloc(4 * textures_size * 2 * sizeof(double));
            memcpy(new_arr, textures, texture_count * 4 * sizeof(double));
            free(textures);
            textures = new_arr;
            textures_itr = new_arr + texture_count * 4;
            textures_size *= 2;
        }
        if (normals_size <= (normal_count + 4)) {
            double *new_arr = (double *)malloc(4 * normals_size * 2 * sizeof(double));
            memcpy(new_arr, normals, normal_count * 4 * sizeof(double));
            free(normals);
            normals = new_arr;
            normals_itr = new_arr + normal_count * 4;
            normals_size *= 2;
        }
        if (groups_size == groups_count) {
            group_with_name *new_groups = (group_with_name *)malloc(2 * groups_size * sizeof(group_with_name));
            memcpy(new_groups, groups, groups_count * sizeof(group_with_name));
            free(groups);
            groups = new_groups;
            groups_size = 2 * groups_size;
        }
    }

    if (groups_count > 1) {
        Shape all_groups = array_of_shapes(groups_count + 1);
        int i;
        for (i = 1; i < groups_count; i++) {
            memcpy(all_groups + i, (groups + i), sizeof(struct shape));
        }
        
        group(all_groups, all_groups + 1, groups_count);

        /* TODO free everthing which is no longer used */

        return all_groups;
    }

    /* TODO free everthing which is no longer used */

    return groups->group;
}
