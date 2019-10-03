#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>

#include "src/libs/canvas/canvas.h"
#include "src/libs/linalg/linalg.h"
#include "src/libs/obj_loader/obj_loader.h"
#include "src/renderer/renderer.h"
#include "src/renderer/camera.h"
#include "src/renderer/world.h"
#include "src/pattern/pattern.h"
#include "src/shapes/shapes.h"
#include "src/shapes/cone.h"
#include "src/shapes/csg.h"
#include "src/shapes/cube.h"
#include "src/shapes/cylinder.h"
#include "src/shapes/group.h"
#include "src/shapes/plane.h"
#include "src/shapes/sphere.h"
#include "src/shapes/triangle.h"
#include "src/shapes/toroid.h"

int main()
{
    /* config */
    size_t thread_count = 4;
//    size_t timeout = 30;
    size_t divide_threshold = 1;
//    bool clamping = false;
    /* end config */

    /* camera */
    struct aperture ap;
    aperture(POINT_APERTURE, 0.0, false, &ap);

    Point from = { 0.0000000000, 2.5000000000, -10.0000000000, 1.0 };
    Point to = { 0.0000000000, 1.0000000000, 0.0000000000, 1.0 };
    Vector up = { 0.0000000000, 1.0000000000, 0.0000000000, 0.0 };
    Matrix camera_xform;
    view_transform(from, to, up, camera_xform);

    Camera cam = camera(1200, 480, 1.2000000000/*field_of_view*/, 1.0000000000/*distance*/, &ap, 1/*sample_num*/, camera_xform);

    /* end camera */

    /* lights */
    Light all_lights = array_of_lights(4);

    /* point light 0 */
    Light point_light_0 = all_lights + 0;
    Point point_light_0_point = { -10.0000000000, 100.0000000000, -100.0000000000, 1.0 };
    Color point_light_0_intensity = color(0.8000000000, 0.8000000000, 0.8000000000);
    point_light(point_light_0_point, point_light_0_intensity, point_light_0);

    /* end point light 0 */
    /* point light 1 */
    Light point_light_1 = all_lights + 1;
    Point point_light_1_point = { 0.0000000000, 100.0000000000, 0.0000000000, 1.0 };
    Color point_light_1_intensity = color(0.1000000000, 0.1000000000, 0.1000000000);
    point_light(point_light_1_point, point_light_1_intensity, point_light_1);

    /* end point light 1 */
    /* point light 2 */
    Light point_light_2 = all_lights + 2;
    Point point_light_2_point = { 100.0000000000, 10.0000000000, -25.0000000000, 1.0 };
    Color point_light_2_intensity = color(0.2000000000, 0.2000000000, 0.2000000000);
    point_light(point_light_2_point, point_light_2_intensity, point_light_2);

    /* end point light 2 */
    /* point light 3 */
    Light point_light_3 = all_lights + 3;
    Point point_light_3_point = { -100.0000000000, 10.0000000000, -25.0000000000, 1.0 };
    Color point_light_3_intensity = color(0.2000000000, 0.2000000000, 0.2000000000);
    point_light(point_light_3_point, point_light_3_intensity, point_light_3);

    /* end point light 3 */

    /* end lights */

    /* shapes */
    Shape all_shapes = array_of_shapes(6);

    /* shape 0 */

    /* children for 0 */
    Shape shape_0_children = array_of_shapes(2);

    
    Pattern pattern_0_child_0 = NULL;

    Material material_0_child_0 = material_alloc();
    material_0_child_0->color[0] = 0.2000000000;
    material_0_child_0->color[1] = 0.2000000000;
    material_0_child_0->color[2] = 0.2000000000;
    material_0_child_0->ambient = 0.0000000000;
    material_0_child_0->diffuse = 0.8000000000;
    material_0_child_0->specular = 0.0000000000;
    material_0_child_0->shininess = 200.0000000000;
    material_0_child_0->reflective = 0.2000000000;
    material_0_child_0->transparency = 0.0000000000;
    material_0_child_0->refractive_index = 1.0000000000;
    material_0_child_0->casts_shadow = true;
    material_0_child_0->pattern = pattern_0_child_0;

    Matrix transform_0_child_0;
    matrix_identity(transform_0_child_0);
    Shape shape_0_child_0 = shape_0_children + 0;
    cylinder(shape_0_child_0);
    shape_set_material(shape_0_child_0, material_0_child_0);
    shape_set_transform(shape_0_child_0, transform_0_child_0);
    shape_0_child_0->fields.cylinder.minimum = -0.1500000000;
    shape_0_child_0->fields.cylinder.maximum = 0.0000000000;
    shape_0_child_0->fields.cylinder.closed = true;


    
    /* children for 0_child_1 */
    Shape shape_0_child_1_children = array_of_shapes(2);

    Pattern pattern_0_child_1_child_0 = NULL;

    Material material_0_child_1_child_0 = material_alloc();
    material_0_child_1_child_0->color[0] = 1.0000000000;
    material_0_child_1_child_0->color[1] = 0.0000000000;
    material_0_child_1_child_0->color[2] = 0.1000000000;
    material_0_child_1_child_0->ambient = 0.1000000000;
    material_0_child_1_child_0->diffuse = 0.6000000000;
    material_0_child_1_child_0->specular = 0.3000000000;
    material_0_child_1_child_0->shininess = 15.0000000000;
    material_0_child_1_child_0->reflective = 0.0000000000;
    material_0_child_1_child_0->transparency = 0.0000000000;
    material_0_child_1_child_0->refractive_index = 1.0000000000;
    material_0_child_1_child_0->casts_shadow = true;
    material_0_child_1_child_0->pattern = pattern_0_child_1_child_0;

    Matrix transform_0_child_1_child_0, transform_0_child_1_child_0_tmp;
    matrix_identity(transform_0_child_1_child_0);
    matrix_translate(0.0000000000, 0.1217000000, 0.0000000000, transform_0_child_1_child_0_tmp);
    transform_chain(transform_0_child_1_child_0_tmp, transform_0_child_1_child_0);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_0_child_1_child_0_tmp);
    transform_chain(transform_0_child_1_child_0_tmp, transform_0_child_1_child_0);

    Shape shape_0_child_1_child_0 = shape_0_child_1_children + 0;
    
    if (access("scenes/bounding_boxes/dragon.obj", F_OK ) == -1 ) {
        printf("file 'scenes/bounding_boxes/dragon.obj' does not exist.");
        return 1;
    }
    printf("Loading resource 'scenes/bounding_boxes/dragon.obj'... ");
    fflush(stdout);
    construct_group_from_obj_file("scenes/bounding_boxes/dragon.obj", shape_0_child_1_child_0);
    printf("Done!\n");
    fflush(stdout);
    shape_set_material_recursive(shape_0_child_1_child_0, material_0_child_1_child_0);
    shape_set_transform(shape_0_child_1_child_0, transform_0_child_1_child_0);


    
    Pattern pattern_0_child_1_child_1 = NULL;

    Material material_0_child_1_child_1 = material_alloc();
    material_0_child_1_child_1->color[0] = 1.0000000000;
    material_0_child_1_child_1->color[1] = 1.0000000000;
    material_0_child_1_child_1->color[2] = 1.0000000000;
    material_0_child_1_child_1->ambient = 0.0000000000;
    material_0_child_1_child_1->diffuse = 0.4000000000;
    material_0_child_1_child_1->specular = 0.0000000000;
    material_0_child_1_child_1->shininess = 200.0000000000;
    material_0_child_1_child_1->reflective = 0.0000000000;
    material_0_child_1_child_1->transparency = 0.6000000000;
    material_0_child_1_child_1->refractive_index = 1.0000000000;
    material_0_child_1_child_1->casts_shadow = false;
    material_0_child_1_child_1->pattern = pattern_0_child_1_child_1;

    Matrix transform_0_child_1_child_1, transform_0_child_1_child_1_tmp;
    matrix_identity(transform_0_child_1_child_1);
    matrix_translate(1.0000000000, 1.0000000000, 1.0000000000, transform_0_child_1_child_1_tmp);
    transform_chain(transform_0_child_1_child_1_tmp, transform_0_child_1_child_1);
    matrix_scale(3.7333500000, 2.5845000000, 1.6283000000, transform_0_child_1_child_1_tmp);
    transform_chain(transform_0_child_1_child_1_tmp, transform_0_child_1_child_1);
    matrix_translate(-3.9863000000, -0.1217000000, -1.1820000000, transform_0_child_1_child_1_tmp);
    transform_chain(transform_0_child_1_child_1_tmp, transform_0_child_1_child_1);
    matrix_translate(0.0000000000, 0.1216900000, 0.0000000000, transform_0_child_1_child_1_tmp);
    transform_chain(transform_0_child_1_child_1_tmp, transform_0_child_1_child_1);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_0_child_1_child_1_tmp);
    transform_chain(transform_0_child_1_child_1_tmp, transform_0_child_1_child_1);

    Shape shape_0_child_1_child_1 = shape_0_child_1_children + 1;
    cube(shape_0_child_1_child_1);
    shape_set_material(shape_0_child_1_child_1, material_0_child_1_child_1);
    shape_set_transform(shape_0_child_1_child_1, transform_0_child_1_child_1);

    /* end children for 0_child_1 */

    Matrix transform_0_child_1;
    matrix_identity(transform_0_child_1);
    Shape shape_0_child_1 = shape_0_children + 1;
    group(shape_0_child_1, shape_0_child_1_children, 2);
    shape_free(shape_0_child_1_children);
    shape_set_transform(shape_0_child_1, transform_0_child_1);

    /* end children for 0 */

    Matrix transform_0;
    matrix_translate(0.0000000000, 2.0000000000, 0.0000000000, transform_0);
    Shape shape_0 = all_shapes + 0;
    group(shape_0, shape_0_children, 2);
    shape_free(shape_0_children);
    shape_set_transform(shape_0, transform_0);

    /* end shape 0 */
    /* shape 1 */

    /* children for 1 */
    Shape shape_1_children = array_of_shapes(2);

    
    Pattern pattern_1_child_0 = NULL;

    Material material_1_child_0 = material_alloc();
    material_1_child_0->color[0] = 0.2000000000;
    material_1_child_0->color[1] = 0.2000000000;
    material_1_child_0->color[2] = 0.2000000000;
    material_1_child_0->ambient = 0.0000000000;
    material_1_child_0->diffuse = 0.8000000000;
    material_1_child_0->specular = 0.0000000000;
    material_1_child_0->shininess = 200.0000000000;
    material_1_child_0->reflective = 0.2000000000;
    material_1_child_0->transparency = 0.0000000000;
    material_1_child_0->refractive_index = 1.0000000000;
    material_1_child_0->casts_shadow = true;
    material_1_child_0->pattern = pattern_1_child_0;

    Matrix transform_1_child_0;
    matrix_identity(transform_1_child_0);
    Shape shape_1_child_0 = shape_1_children + 0;
    cylinder(shape_1_child_0);
    shape_set_material(shape_1_child_0, material_1_child_0);
    shape_set_transform(shape_1_child_0, transform_1_child_0);
    shape_1_child_0->fields.cylinder.minimum = -0.1500000000;
    shape_1_child_0->fields.cylinder.maximum = 0.0000000000;
    shape_1_child_0->fields.cylinder.closed = true;


    
    /* children for 1_child_1 */
    Shape shape_1_child_1_children = array_of_shapes(2);

    Pattern pattern_1_child_1_child_0 = NULL;

    Material material_1_child_1_child_0 = material_alloc();
    material_1_child_1_child_0->color[0] = 1.0000000000;
    material_1_child_1_child_0->color[1] = 0.5000000000;
    material_1_child_1_child_0->color[2] = 0.1000000000;
    material_1_child_1_child_0->ambient = 0.1000000000;
    material_1_child_1_child_0->diffuse = 0.6000000000;
    material_1_child_1_child_0->specular = 0.3000000000;
    material_1_child_1_child_0->shininess = 15.0000000000;
    material_1_child_1_child_0->reflective = 0.0000000000;
    material_1_child_1_child_0->transparency = 0.0000000000;
    material_1_child_1_child_0->refractive_index = 1.0000000000;
    material_1_child_1_child_0->casts_shadow = true;
    material_1_child_1_child_0->pattern = pattern_1_child_1_child_0;

    Matrix transform_1_child_1_child_0, transform_1_child_1_child_0_tmp;
    matrix_identity(transform_1_child_1_child_0);
    matrix_translate(0.0000000000, 0.1217000000, 0.0000000000, transform_1_child_1_child_0_tmp);
    transform_chain(transform_1_child_1_child_0_tmp, transform_1_child_1_child_0);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_1_child_1_child_0_tmp);
    transform_chain(transform_1_child_1_child_0_tmp, transform_1_child_1_child_0);

    Shape shape_1_child_1_child_0 = shape_1_child_1_children + 0;
    
    if (access("scenes/bounding_boxes/dragon.obj", F_OK ) == -1 ) {
        printf("file 'scenes/bounding_boxes/dragon.obj' does not exist.");
        return 1;
    }
    printf("Loading resource 'scenes/bounding_boxes/dragon.obj'... ");
    fflush(stdout);
    construct_group_from_obj_file("scenes/bounding_boxes/dragon.obj", shape_1_child_1_child_0);
    printf("Done!\n");
    fflush(stdout);
    shape_set_material_recursive(shape_1_child_1_child_0, material_1_child_1_child_0);
    shape_set_transform(shape_1_child_1_child_0, transform_1_child_1_child_0);


    
    Pattern pattern_1_child_1_child_1 = NULL;

    Material material_1_child_1_child_1 = material_alloc();
    material_1_child_1_child_1->color[0] = 1.0000000000;
    material_1_child_1_child_1->color[1] = 1.0000000000;
    material_1_child_1_child_1->color[2] = 1.0000000000;
    material_1_child_1_child_1->ambient = 0.0000000000;
    material_1_child_1_child_1->diffuse = 0.2000000000;
    material_1_child_1_child_1->specular = 0.0000000000;
    material_1_child_1_child_1->shininess = 200.0000000000;
    material_1_child_1_child_1->reflective = 0.0000000000;
    material_1_child_1_child_1->transparency = 0.8000000000;
    material_1_child_1_child_1->refractive_index = 1.0000000000;
    material_1_child_1_child_1->casts_shadow = false;
    material_1_child_1_child_1->pattern = pattern_1_child_1_child_1;

    Matrix transform_1_child_1_child_1, transform_1_child_1_child_1_tmp;
    matrix_identity(transform_1_child_1_child_1);
    matrix_translate(1.0000000000, 1.0000000000, 1.0000000000, transform_1_child_1_child_1_tmp);
    transform_chain(transform_1_child_1_child_1_tmp, transform_1_child_1_child_1);
    matrix_scale(3.7333500000, 2.5845000000, 1.6283000000, transform_1_child_1_child_1_tmp);
    transform_chain(transform_1_child_1_child_1_tmp, transform_1_child_1_child_1);
    matrix_translate(-3.9863000000, -0.1217000000, -1.1820000000, transform_1_child_1_child_1_tmp);
    transform_chain(transform_1_child_1_child_1_tmp, transform_1_child_1_child_1);
    matrix_translate(0.0000000000, 0.1216900000, 0.0000000000, transform_1_child_1_child_1_tmp);
    transform_chain(transform_1_child_1_child_1_tmp, transform_1_child_1_child_1);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_1_child_1_child_1_tmp);
    transform_chain(transform_1_child_1_child_1_tmp, transform_1_child_1_child_1);

    Shape shape_1_child_1_child_1 = shape_1_child_1_children + 1;
    cube(shape_1_child_1_child_1);
    shape_set_material(shape_1_child_1_child_1, material_1_child_1_child_1);
    shape_set_transform(shape_1_child_1_child_1, transform_1_child_1_child_1);

    /* end children for 1_child_1 */

    Matrix transform_1_child_1, transform_1_child_1_tmp;
    matrix_identity(transform_1_child_1);
    matrix_rotate_y(4.0000000000, transform_1_child_1_tmp);
    transform_chain(transform_1_child_1_tmp, transform_1_child_1);
    matrix_scale(0.7500000000, 0.7500000000, 0.7500000000, transform_1_child_1_tmp);
    transform_chain(transform_1_child_1_tmp, transform_1_child_1);

    Shape shape_1_child_1 = shape_1_children + 1;
    group(shape_1_child_1, shape_1_child_1_children, 2);
    shape_free(shape_1_child_1_children);
    shape_set_transform(shape_1_child_1, transform_1_child_1);

    /* end children for 1 */

    Matrix transform_1;
    matrix_translate(2.0000000000, 1.0000000000, -1.0000000000, transform_1);
    Shape shape_1 = all_shapes + 1;
    group(shape_1, shape_1_children, 2);
    shape_free(shape_1_children);
    shape_set_transform(shape_1, transform_1);

    /* end shape 1 */
    /* shape 2 */

    /* children for 2 */
    Shape shape_2_children = array_of_shapes(2);

    
    Pattern pattern_2_child_0 = NULL;

    Material material_2_child_0 = material_alloc();
    material_2_child_0->color[0] = 0.2000000000;
    material_2_child_0->color[1] = 0.2000000000;
    material_2_child_0->color[2] = 0.2000000000;
    material_2_child_0->ambient = 0.0000000000;
    material_2_child_0->diffuse = 0.8000000000;
    material_2_child_0->specular = 0.0000000000;
    material_2_child_0->shininess = 200.0000000000;
    material_2_child_0->reflective = 0.2000000000;
    material_2_child_0->transparency = 0.0000000000;
    material_2_child_0->refractive_index = 1.0000000000;
    material_2_child_0->casts_shadow = true;
    material_2_child_0->pattern = pattern_2_child_0;

    Matrix transform_2_child_0;
    matrix_identity(transform_2_child_0);
    Shape shape_2_child_0 = shape_2_children + 0;
    cylinder(shape_2_child_0);
    shape_set_material(shape_2_child_0, material_2_child_0);
    shape_set_transform(shape_2_child_0, transform_2_child_0);
    shape_2_child_0->fields.cylinder.minimum = -0.1500000000;
    shape_2_child_0->fields.cylinder.maximum = 0.0000000000;
    shape_2_child_0->fields.cylinder.closed = true;


    
    /* children for 2_child_1 */
    Shape shape_2_child_1_children = array_of_shapes(2);

    Pattern pattern_2_child_1_child_0 = NULL;

    Material material_2_child_1_child_0 = material_alloc();
    material_2_child_1_child_0->color[0] = 0.9000000000;
    material_2_child_1_child_0->color[1] = 0.5000000000;
    material_2_child_1_child_0->color[2] = 0.1000000000;
    material_2_child_1_child_0->ambient = 0.1000000000;
    material_2_child_1_child_0->diffuse = 0.6000000000;
    material_2_child_1_child_0->specular = 0.3000000000;
    material_2_child_1_child_0->shininess = 15.0000000000;
    material_2_child_1_child_0->reflective = 0.0000000000;
    material_2_child_1_child_0->transparency = 0.0000000000;
    material_2_child_1_child_0->refractive_index = 1.0000000000;
    material_2_child_1_child_0->casts_shadow = true;
    material_2_child_1_child_0->pattern = pattern_2_child_1_child_0;

    Matrix transform_2_child_1_child_0, transform_2_child_1_child_0_tmp;
    matrix_identity(transform_2_child_1_child_0);
    matrix_translate(0.0000000000, 0.1217000000, 0.0000000000, transform_2_child_1_child_0_tmp);
    transform_chain(transform_2_child_1_child_0_tmp, transform_2_child_1_child_0);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_2_child_1_child_0_tmp);
    transform_chain(transform_2_child_1_child_0_tmp, transform_2_child_1_child_0);

    Shape shape_2_child_1_child_0 = shape_2_child_1_children + 0;
    
    if (access("scenes/bounding_boxes/dragon.obj", F_OK ) == -1 ) {
        printf("file 'scenes/bounding_boxes/dragon.obj' does not exist.");
        return 1;
    }
    printf("Loading resource 'scenes/bounding_boxes/dragon.obj'... ");
    fflush(stdout);
    construct_group_from_obj_file("scenes/bounding_boxes/dragon.obj", shape_2_child_1_child_0);
    printf("Done!\n");
    fflush(stdout);
    shape_set_material_recursive(shape_2_child_1_child_0, material_2_child_1_child_0);
    shape_set_transform(shape_2_child_1_child_0, transform_2_child_1_child_0);


    
    Pattern pattern_2_child_1_child_1 = NULL;

    Material material_2_child_1_child_1 = material_alloc();
    material_2_child_1_child_1->color[0] = 1.0000000000;
    material_2_child_1_child_1->color[1] = 1.0000000000;
    material_2_child_1_child_1->color[2] = 1.0000000000;
    material_2_child_1_child_1->ambient = 0.0000000000;
    material_2_child_1_child_1->diffuse = 0.2000000000;
    material_2_child_1_child_1->specular = 0.0000000000;
    material_2_child_1_child_1->shininess = 200.0000000000;
    material_2_child_1_child_1->reflective = 0.0000000000;
    material_2_child_1_child_1->transparency = 0.8000000000;
    material_2_child_1_child_1->refractive_index = 1.0000000000;
    material_2_child_1_child_1->casts_shadow = false;
    material_2_child_1_child_1->pattern = pattern_2_child_1_child_1;

    Matrix transform_2_child_1_child_1, transform_2_child_1_child_1_tmp;
    matrix_identity(transform_2_child_1_child_1);
    matrix_translate(1.0000000000, 1.0000000000, 1.0000000000, transform_2_child_1_child_1_tmp);
    transform_chain(transform_2_child_1_child_1_tmp, transform_2_child_1_child_1);
    matrix_scale(3.7333500000, 2.5845000000, 1.6283000000, transform_2_child_1_child_1_tmp);
    transform_chain(transform_2_child_1_child_1_tmp, transform_2_child_1_child_1);
    matrix_translate(-3.9863000000, -0.1217000000, -1.1820000000, transform_2_child_1_child_1_tmp);
    transform_chain(transform_2_child_1_child_1_tmp, transform_2_child_1_child_1);
    matrix_translate(0.0000000000, 0.1216900000, 0.0000000000, transform_2_child_1_child_1_tmp);
    transform_chain(transform_2_child_1_child_1_tmp, transform_2_child_1_child_1);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_2_child_1_child_1_tmp);
    transform_chain(transform_2_child_1_child_1_tmp, transform_2_child_1_child_1);

    Shape shape_2_child_1_child_1 = shape_2_child_1_children + 1;
    cube(shape_2_child_1_child_1);
    shape_set_material(shape_2_child_1_child_1, material_2_child_1_child_1);
    shape_set_transform(shape_2_child_1_child_1, transform_2_child_1_child_1);

    /* end children for 2_child_1 */

    Matrix transform_2_child_1, transform_2_child_1_tmp;
    matrix_identity(transform_2_child_1);
    matrix_rotate_y(-0.4000000000, transform_2_child_1_tmp);
    transform_chain(transform_2_child_1_tmp, transform_2_child_1);
    matrix_scale(0.7500000000, 0.7500000000, 0.7500000000, transform_2_child_1_tmp);
    transform_chain(transform_2_child_1_tmp, transform_2_child_1);

    Shape shape_2_child_1 = shape_2_children + 1;
    group(shape_2_child_1, shape_2_child_1_children, 2);
    shape_free(shape_2_child_1_children);
    shape_set_transform(shape_2_child_1, transform_2_child_1);

    /* end children for 2 */

    Matrix transform_2;
    matrix_translate(-2.0000000000, 0.7500000000, -1.0000000000, transform_2);
    Shape shape_2 = all_shapes + 2;
    group(shape_2, shape_2_children, 2);
    shape_free(shape_2_children);
    shape_set_transform(shape_2, transform_2);

    /* end shape 2 */
    /* shape 3 */

    /* children for 3 */
    Shape shape_3_children = array_of_shapes(2);

    
    Pattern pattern_3_child_0 = NULL;

    Material material_3_child_0 = material_alloc();
    material_3_child_0->color[0] = 0.2000000000;
    material_3_child_0->color[1] = 0.2000000000;
    material_3_child_0->color[2] = 0.2000000000;
    material_3_child_0->ambient = 0.0000000000;
    material_3_child_0->diffuse = 0.8000000000;
    material_3_child_0->specular = 0.0000000000;
    material_3_child_0->shininess = 200.0000000000;
    material_3_child_0->reflective = 0.2000000000;
    material_3_child_0->transparency = 0.0000000000;
    material_3_child_0->refractive_index = 1.0000000000;
    material_3_child_0->casts_shadow = true;
    material_3_child_0->pattern = pattern_3_child_0;

    Matrix transform_3_child_0;
    matrix_identity(transform_3_child_0);
    Shape shape_3_child_0 = shape_3_children + 0;
    cylinder(shape_3_child_0);
    shape_set_material(shape_3_child_0, material_3_child_0);
    shape_set_transform(shape_3_child_0, transform_3_child_0);
    shape_3_child_0->fields.cylinder.minimum = -0.1500000000;
    shape_3_child_0->fields.cylinder.maximum = 0.0000000000;
    shape_3_child_0->fields.cylinder.closed = true;


    
    /* children for 3_child_1 */
    Shape shape_3_child_1_children = array_of_shapes(2);

    Pattern pattern_3_child_1_child_0 = NULL;

    Material material_3_child_1_child_0 = material_alloc();
    material_3_child_1_child_0->color[0] = 1.0000000000;
    material_3_child_1_child_0->color[1] = 0.9000000000;
    material_3_child_1_child_0->color[2] = 0.1000000000;
    material_3_child_1_child_0->ambient = 0.1000000000;
    material_3_child_1_child_0->diffuse = 0.6000000000;
    material_3_child_1_child_0->specular = 0.3000000000;
    material_3_child_1_child_0->shininess = 15.0000000000;
    material_3_child_1_child_0->reflective = 0.0000000000;
    material_3_child_1_child_0->transparency = 0.0000000000;
    material_3_child_1_child_0->refractive_index = 1.0000000000;
    material_3_child_1_child_0->casts_shadow = true;
    material_3_child_1_child_0->pattern = pattern_3_child_1_child_0;

    Matrix transform_3_child_1_child_0, transform_3_child_1_child_0_tmp;
    matrix_identity(transform_3_child_1_child_0);
    matrix_translate(0.0000000000, 0.1217000000, 0.0000000000, transform_3_child_1_child_0_tmp);
    transform_chain(transform_3_child_1_child_0_tmp, transform_3_child_1_child_0);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_3_child_1_child_0_tmp);
    transform_chain(transform_3_child_1_child_0_tmp, transform_3_child_1_child_0);

    Shape shape_3_child_1_child_0 = shape_3_child_1_children + 0;
    
    if (access("scenes/bounding_boxes/dragon.obj", F_OK ) == -1 ) {
        printf("file 'scenes/bounding_boxes/dragon.obj' does not exist.");
        return 1;
    }
    printf("Loading resource 'scenes/bounding_boxes/dragon.obj'... ");
    fflush(stdout);
    construct_group_from_obj_file("scenes/bounding_boxes/dragon.obj", shape_3_child_1_child_0);
    printf("Done!\n");
    fflush(stdout);
    shape_set_material_recursive(shape_3_child_1_child_0, material_3_child_1_child_0);
    shape_set_transform(shape_3_child_1_child_0, transform_3_child_1_child_0);


    
    Pattern pattern_3_child_1_child_1 = NULL;

    Material material_3_child_1_child_1 = material_alloc();
    material_3_child_1_child_1->color[0] = 1.0000000000;
    material_3_child_1_child_1->color[1] = 1.0000000000;
    material_3_child_1_child_1->color[2] = 1.0000000000;
    material_3_child_1_child_1->ambient = 0.0000000000;
    material_3_child_1_child_1->diffuse = 0.1000000000;
    material_3_child_1_child_1->specular = 0.0000000000;
    material_3_child_1_child_1->shininess = 200.0000000000;
    material_3_child_1_child_1->reflective = 0.0000000000;
    material_3_child_1_child_1->transparency = 0.9000000000;
    material_3_child_1_child_1->refractive_index = 1.0000000000;
    material_3_child_1_child_1->casts_shadow = false;
    material_3_child_1_child_1->pattern = pattern_3_child_1_child_1;

    Matrix transform_3_child_1_child_1, transform_3_child_1_child_1_tmp;
    matrix_identity(transform_3_child_1_child_1);
    matrix_translate(1.0000000000, 1.0000000000, 1.0000000000, transform_3_child_1_child_1_tmp);
    transform_chain(transform_3_child_1_child_1_tmp, transform_3_child_1_child_1);
    matrix_scale(3.7333500000, 2.5845000000, 1.6283000000, transform_3_child_1_child_1_tmp);
    transform_chain(transform_3_child_1_child_1_tmp, transform_3_child_1_child_1);
    matrix_translate(-3.9863000000, -0.1217000000, -1.1820000000, transform_3_child_1_child_1_tmp);
    transform_chain(transform_3_child_1_child_1_tmp, transform_3_child_1_child_1);
    matrix_translate(0.0000000000, 0.1216900000, 0.0000000000, transform_3_child_1_child_1_tmp);
    transform_chain(transform_3_child_1_child_1_tmp, transform_3_child_1_child_1);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_3_child_1_child_1_tmp);
    transform_chain(transform_3_child_1_child_1_tmp, transform_3_child_1_child_1);

    Shape shape_3_child_1_child_1 = shape_3_child_1_children + 1;
    cube(shape_3_child_1_child_1);
    shape_set_material(shape_3_child_1_child_1, material_3_child_1_child_1);
    shape_set_transform(shape_3_child_1_child_1, transform_3_child_1_child_1);

    /* end children for 3_child_1 */

    Matrix transform_3_child_1, transform_3_child_1_tmp;
    matrix_identity(transform_3_child_1);
    matrix_rotate_y(-0.2000000000, transform_3_child_1_tmp);
    transform_chain(transform_3_child_1_tmp, transform_3_child_1);
    matrix_scale(0.5000000000, 0.5000000000, 0.5000000000, transform_3_child_1_tmp);
    transform_chain(transform_3_child_1_tmp, transform_3_child_1);

    Shape shape_3_child_1 = shape_3_children + 1;
    group(shape_3_child_1, shape_3_child_1_children, 2);
    shape_free(shape_3_child_1_children);
    shape_set_transform(shape_3_child_1, transform_3_child_1);

    /* end children for 3 */

    Matrix transform_3;
    matrix_translate(-4.0000000000, 0.0000000000, -2.0000000000, transform_3);
    Shape shape_3 = all_shapes + 3;
    group(shape_3, shape_3_children, 2);
    shape_free(shape_3_children);
    shape_set_transform(shape_3, transform_3);

    /* end shape 3 */
    /* shape 4 */

    /* children for 4 */
    Shape shape_4_children = array_of_shapes(2);

    
    Pattern pattern_4_child_0 = NULL;

    Material material_4_child_0 = material_alloc();
    material_4_child_0->color[0] = 0.2000000000;
    material_4_child_0->color[1] = 0.2000000000;
    material_4_child_0->color[2] = 0.2000000000;
    material_4_child_0->ambient = 0.0000000000;
    material_4_child_0->diffuse = 0.8000000000;
    material_4_child_0->specular = 0.0000000000;
    material_4_child_0->shininess = 200.0000000000;
    material_4_child_0->reflective = 0.2000000000;
    material_4_child_0->transparency = 0.0000000000;
    material_4_child_0->refractive_index = 1.0000000000;
    material_4_child_0->casts_shadow = true;
    material_4_child_0->pattern = pattern_4_child_0;

    Matrix transform_4_child_0;
    matrix_identity(transform_4_child_0);
    Shape shape_4_child_0 = shape_4_children + 0;
    cylinder(shape_4_child_0);
    shape_set_material(shape_4_child_0, material_4_child_0);
    shape_set_transform(shape_4_child_0, transform_4_child_0);
    shape_4_child_0->fields.cylinder.minimum = -0.1500000000;
    shape_4_child_0->fields.cylinder.maximum = 0.0000000000;
    shape_4_child_0->fields.cylinder.closed = true;


    
    /* children for 4_child_1 */
    Shape shape_4_child_1_children = array_of_shapes(2);

    Pattern pattern_4_child_1_child_0 = NULL;

    Material material_4_child_1_child_0 = material_alloc();
    material_4_child_1_child_0->color[0] = 0.9000000000;
    material_4_child_1_child_0->color[1] = 1.0000000000;
    material_4_child_1_child_0->color[2] = 0.1000000000;
    material_4_child_1_child_0->ambient = 0.1000000000;
    material_4_child_1_child_0->diffuse = 0.6000000000;
    material_4_child_1_child_0->specular = 0.3000000000;
    material_4_child_1_child_0->shininess = 15.0000000000;
    material_4_child_1_child_0->reflective = 0.0000000000;
    material_4_child_1_child_0->transparency = 0.0000000000;
    material_4_child_1_child_0->refractive_index = 1.0000000000;
    material_4_child_1_child_0->casts_shadow = true;
    material_4_child_1_child_0->pattern = pattern_4_child_1_child_0;

    Matrix transform_4_child_1_child_0, transform_4_child_1_child_0_tmp;
    matrix_identity(transform_4_child_1_child_0);
    matrix_translate(0.0000000000, 0.1217000000, 0.0000000000, transform_4_child_1_child_0_tmp);
    transform_chain(transform_4_child_1_child_0_tmp, transform_4_child_1_child_0);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_4_child_1_child_0_tmp);
    transform_chain(transform_4_child_1_child_0_tmp, transform_4_child_1_child_0);

    Shape shape_4_child_1_child_0 = shape_4_child_1_children + 0;
    
    if (access("scenes/bounding_boxes/dragon.obj", F_OK ) == -1 ) {
        printf("file 'scenes/bounding_boxes/dragon.obj' does not exist.");
        return 1;
    }
    printf("Loading resource 'scenes/bounding_boxes/dragon.obj'... ");
    fflush(stdout);
    construct_group_from_obj_file("scenes/bounding_boxes/dragon.obj", shape_4_child_1_child_0);
    printf("Done!\n");
    fflush(stdout);
    shape_set_material_recursive(shape_4_child_1_child_0, material_4_child_1_child_0);
    shape_set_transform(shape_4_child_1_child_0, transform_4_child_1_child_0);


    
    Pattern pattern_4_child_1_child_1 = NULL;

    Material material_4_child_1_child_1 = material_alloc();
    material_4_child_1_child_1->color[0] = 1.0000000000;
    material_4_child_1_child_1->color[1] = 1.0000000000;
    material_4_child_1_child_1->color[2] = 1.0000000000;
    material_4_child_1_child_1->ambient = 0.0000000000;
    material_4_child_1_child_1->diffuse = 0.1000000000;
    material_4_child_1_child_1->specular = 0.0000000000;
    material_4_child_1_child_1->shininess = 200.0000000000;
    material_4_child_1_child_1->reflective = 0.0000000000;
    material_4_child_1_child_1->transparency = 0.9000000000;
    material_4_child_1_child_1->refractive_index = 1.0000000000;
    material_4_child_1_child_1->casts_shadow = false;
    material_4_child_1_child_1->pattern = pattern_4_child_1_child_1;

    Matrix transform_4_child_1_child_1, transform_4_child_1_child_1_tmp;
    matrix_identity(transform_4_child_1_child_1);
    matrix_translate(1.0000000000, 1.0000000000, 1.0000000000, transform_4_child_1_child_1_tmp);
    transform_chain(transform_4_child_1_child_1_tmp, transform_4_child_1_child_1);
    matrix_scale(3.7333500000, 2.5845000000, 1.6283000000, transform_4_child_1_child_1_tmp);
    transform_chain(transform_4_child_1_child_1_tmp, transform_4_child_1_child_1);
    matrix_translate(-3.9863000000, -0.1217000000, -1.1820000000, transform_4_child_1_child_1_tmp);
    transform_chain(transform_4_child_1_child_1_tmp, transform_4_child_1_child_1);
    matrix_translate(0.0000000000, 0.1216900000, 0.0000000000, transform_4_child_1_child_1_tmp);
    transform_chain(transform_4_child_1_child_1_tmp, transform_4_child_1_child_1);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_4_child_1_child_1_tmp);
    transform_chain(transform_4_child_1_child_1_tmp, transform_4_child_1_child_1);

    Shape shape_4_child_1_child_1 = shape_4_child_1_children + 1;
    cube(shape_4_child_1_child_1);
    shape_set_material(shape_4_child_1_child_1, material_4_child_1_child_1);
    shape_set_transform(shape_4_child_1_child_1, transform_4_child_1_child_1);

    /* end children for 4_child_1 */

    Matrix transform_4_child_1, transform_4_child_1_tmp;
    matrix_identity(transform_4_child_1);
    matrix_rotate_y(3.3000000000, transform_4_child_1_tmp);
    transform_chain(transform_4_child_1_tmp, transform_4_child_1);
    matrix_scale(0.5000000000, 0.5000000000, 0.5000000000, transform_4_child_1_tmp);
    transform_chain(transform_4_child_1_tmp, transform_4_child_1);

    Shape shape_4_child_1 = shape_4_children + 1;
    group(shape_4_child_1, shape_4_child_1_children, 2);
    shape_free(shape_4_child_1_children);
    shape_set_transform(shape_4_child_1, transform_4_child_1);

    /* end children for 4 */

    Matrix transform_4;
    matrix_translate(4.0000000000, 0.0000000000, -2.0000000000, transform_4);
    Shape shape_4 = all_shapes + 4;
    group(shape_4, shape_4_children, 2);
    shape_free(shape_4_children);
    shape_set_transform(shape_4, transform_4);

    /* end shape 4 */
    /* shape 5 */

    /* children for 5 */
    Shape shape_5_children = array_of_shapes(2);

    
    Pattern pattern_5_child_0 = NULL;

    Material material_5_child_0 = material_alloc();
    material_5_child_0->color[0] = 0.2000000000;
    material_5_child_0->color[1] = 0.2000000000;
    material_5_child_0->color[2] = 0.2000000000;
    material_5_child_0->ambient = 0.0000000000;
    material_5_child_0->diffuse = 0.8000000000;
    material_5_child_0->specular = 0.0000000000;
    material_5_child_0->shininess = 200.0000000000;
    material_5_child_0->reflective = 0.2000000000;
    material_5_child_0->transparency = 0.0000000000;
    material_5_child_0->refractive_index = 1.0000000000;
    material_5_child_0->casts_shadow = true;
    material_5_child_0->pattern = pattern_5_child_0;

    Matrix transform_5_child_0;
    matrix_identity(transform_5_child_0);
    Shape shape_5_child_0 = shape_5_children + 0;
    cylinder(shape_5_child_0);
    shape_set_material(shape_5_child_0, material_5_child_0);
    shape_set_transform(shape_5_child_0, transform_5_child_0);
    shape_5_child_0->fields.cylinder.minimum = -0.1500000000;
    shape_5_child_0->fields.cylinder.maximum = 0.0000000000;
    shape_5_child_0->fields.cylinder.closed = true;


    Pattern pattern_5_child_1 = NULL;

    Material material_5_child_1 = material_alloc();
    material_5_child_1->color[0] = 1.0000000000;
    material_5_child_1->color[1] = 1.0000000000;
    material_5_child_1->color[2] = 1.0000000000;
    material_5_child_1->ambient = 0.1000000000;
    material_5_child_1->diffuse = 0.6000000000;
    material_5_child_1->specular = 0.3000000000;
    material_5_child_1->shininess = 15.0000000000;
    material_5_child_1->reflective = 0.0000000000;
    material_5_child_1->transparency = 0.0000000000;
    material_5_child_1->refractive_index = 1.0000000000;
    material_5_child_1->casts_shadow = true;
    material_5_child_1->pattern = pattern_5_child_1;

    Matrix transform_5_child_1, transform_5_child_1_tmp;
    matrix_identity(transform_5_child_1);
    matrix_translate(0.0000000000, 0.1217000000, 0.0000000000, transform_5_child_1_tmp);
    transform_chain(transform_5_child_1_tmp, transform_5_child_1);
    matrix_scale(0.2680000000, 0.2680000000, 0.2680000000, transform_5_child_1_tmp);
    transform_chain(transform_5_child_1_tmp, transform_5_child_1);
    matrix_rotate_y(3.1415000000, transform_5_child_1_tmp);
    transform_chain(transform_5_child_1_tmp, transform_5_child_1);

    Shape shape_5_child_1 = shape_5_children + 1;
    
    if (access("scenes/bounding_boxes/dragon.obj", F_OK ) == -1 ) {
        printf("file 'scenes/bounding_boxes/dragon.obj' does not exist.");
        return 1;
    }
    printf("Loading resource 'scenes/bounding_boxes/dragon.obj'... ");
    fflush(stdout);
    construct_group_from_obj_file("scenes/bounding_boxes/dragon.obj", shape_5_child_1);
    printf("Done!\n");
    fflush(stdout);
    shape_set_material_recursive(shape_5_child_1, material_5_child_1);
    shape_set_transform(shape_5_child_1, transform_5_child_1);

    /* end children for 5 */

    Matrix transform_5;
    matrix_translate(0.0000000000, 0.5000000000, -4.0000000000, transform_5);
    Shape shape_5 = all_shapes + 5;
    group(shape_5, shape_5_children, 2);
    shape_free(shape_5_children);
    shape_set_transform(shape_5, transform_5);

    /* end shape 5 */
    /* end shapes */

    Shape world_group = array_of_shapes(1);
    group(world_group, all_shapes, 6);
    printf("Balancing scene...");
    fflush(stdout);
    world_group->divide(world_group, divide_threshold);
    printf("Done!\n");
    fflush(stdout);

    World w = world();
    w->lights = all_lights;
    w->lights_num = 4;
    w->shapes = world_group;
    w->shapes_num = 1;
    Canvas c = render_multi(cam, w, cam->sample_num/*usteps*/, cam->sample_num/*vsteps*/, cam->aperture.jitter, thread_count);
    Ppm ppm = construct_ppm(c, true);

    FILE * pFile;
    pFile = fopen ("/tmp/unclamped.ppm", "wb");
    fwrite (ppm->arr, sizeof(unsigned char), ppm->len, pFile);
    fclose (pFile);

    ppm_free(ppm);
    ppm = construct_ppm(c, false);

    pFile = fopen("/tmp/clamped.ppm", "wb");
    fwrite (ppm->arr, sizeof(unsigned char), ppm->len, pFile);
    fclose (pFile);

    ppm_free(ppm);
    canvas_free(c);

    return 0;
}

