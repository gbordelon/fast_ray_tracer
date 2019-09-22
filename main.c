#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "linalg.h"
#include "canvas.h"
#include "renderer.h"
#include "shapes.h"

void
matrix_print_helper(Matrix m, const char * name) {
    static char str_buf[256];
    matrix_to_string(str_buf, sizeof(str_buf), m);
    printf("%s: %s\n", name, str_buf);
}

void
point_print_helper(Point p, const char * name) {
    static char str_buf[256];
    point_to_string(str_buf, sizeof(str_buf), p);
    printf("%s: %s\n", name, str_buf);
}

void
vector_print_helper(Vector v, const char * name) {
    static char str_buf[256];
    vector_to_string(str_buf, sizeof(str_buf), v);
    printf("%s: %s\n", name, str_buf);
}


Canvas
generate_test_image()
{
    Canvas c = canvas_alloc(500, 100);
    Color color = color_default();

    int i, j;
    for (i = 0; i < c->width; i += 1) {
        for (j = 0; j < c->height; j += 1) {
            color->arr[0] = (double)i * (double)j/ 50000.0;
            color->arr[2] = 1.0 - color->arr[0];
            printf("%f %f %f\n", color->arr[0], color->arr[1], color->arr[2]);
            canvas_write_pixel(c, i, j, color);
        }
    }

    color_free(color);
    return c;
}

int
main()
{
/*
    Vector v1 = vector(1.0,2.0,3.0);
    Point p1 = point(1.0, 2.2, 0.3);
    Point p2 = point(0.0,0.5,0.5);
    Vector v2 = vector_from_points_alloc(p1, p2);

    vector_print_helper(v1, "v1");
    vector_print_helper(v2, "v2");
    point_print_helper(p1, "p1");

    printf("v1 magnitude: %f\n", vector_magnitude(v1));

    Vector v3 = vector_normalize_alloc(v1);
    vector_print_helper(v3, "v1 normalize alloc");

    printf("v1 dot v2: %f\n", vector_dot(v1, v2));

    Vector v4 = vector_cross_alloc(v1, v2);
    vector_print_helper(v4, "v1 cross v2");

    v2 = vector_default();
    vector_normalize(v1, v2);
    vector_print_helper(v2, "v1 normalize");

    Matrix m0 = matrix_default();

    matrix_identity(m0);
    matrix_print_helper(m0, "m1");

    matrix_translate(1,2,3, m0);
    matrix_print_helper(m0, "m2");

    Matrix m3 = matrix_scale_alloc(1,2,3);
    matrix_print_helper(m3, "m3");

    m0 = matrix_default();
    matrix_rotate_x(M_PI_4, m0);
    matrix_print_helper(m0, "m4");

    matrix_identity(m0);
    matrix_rotate_y(M_PI_4, m0);
    matrix_print_helper(m0, "m5");

    matrix_identity(m0);
    matrix_rotate_z(M_PI_4, m0);
    matrix_print_helper(m0, "m6");

    matrix_identity(m0);
    matrix_shear(2,4,6,8,10,12, m0);
    matrix_print_helper(m0, "m7");

    matrix_identity(m0);
    matrix_rotate_z(M_PI_4, m0);
    matrix_print_helper(m0, "m8");

    Matrix m9 = matrix_transpose_alloc(m0);
    matrix_print_helper(m9, "m9");

    matrix(-5,2,6,-8, 1,-5,1,8, 7,7,-6,-7, 1,-3,7,4, m0);
    matrix_inverse(m0, m9);
    matrix_print_helper(m9, "identity inverse");

    matrix(3,-9,7,3, 3,-8,2,-9, -4,4,4,1, -6,5,-1,1, m0);
    matrix(8,2,2,2, 3,-1,7,0, 7,0,5,4, 6,-2,0,5, m9);
    Matrix m9_inv = matrix_inverse_alloc(m9);
    Matrix prod = matrix_multiply_alloc(m0,m9);
    Matrix prod2 = matrix_multiply_alloc(prod, m9_inv);
    matrix_print_helper(m0, "a");
    matrix_print_helper(m9, "b");
    matrix_print_helper(m9_inv, "b_inv");
    matrix_print_helper(prod, "c");
    matrix_print_helper(prod2, "d");
    printf("compare a,d: %d\n", matrix_equal(m0,prod2));
    
*/
//    Canvas c = construct_canvas_from_ppm_file("/tmp/myfile2.ppm");
//    Canvas c = generate_test_image();

/*    Point from = point(0.8,0.5,-.8);
    Point to = point(0,0.3,0);
    Vector up = vector(0, 1, 0);
    Camera cam = camera(6, 5, M_PI/40.0, view_transform(from, to, up));
*/
    Point from = point(0,0.5,-4.0);
    Point to = point(0.6,0.0,0);
    Vector up = vector(0, 1, 0);
    Camera cam = camera(400, 400, 1, view_transform(from, to, up));

    point_free(from);
    point_free(to);
    vector_free(up);

    World w = default_world();

    Canvas c = render(cam, w);

    Ppm ppm = construct_ppm(c, true);

    FILE * pFile;
    pFile = fopen ("/tmp/left.ppm", "wb");
    fwrite (ppm->arr, sizeof(unsigned char), ppm->len, pFile);
    fclose (pFile);



/*
    from = point(3.9,3.5,-4.0);
    to = point(0,0.0,0);
    up = vector(0, 1, 0);
    cam = camera(800, 800, .3, view_transform(from, to, up));

    c = render(cam, w);




    Ppm ppm2 = construct_ppm(c, true);

    pFile = fopen ("/tmp/right.ppm", "wb");
    fwrite (ppm2->arr, sizeof(unsigned char), ppm2->len, pFile);
    fclose (pFile);
*/
    printf("file written\n");

    return 0;
}
