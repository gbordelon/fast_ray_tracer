#ifndef LINALG
#define LINALG

#include <stdio.h>
#include <stddef.h>

#define EPSILON 0.00001
#define equal(a,b) (fabs((a) - (b)) < EPSILON)

#define point_equal(a,b) equal((a)->arr[0], (b)->arr[0]) && equal((a)->arr[1], (b)->arr[1]) && equal((a)->arr[2], (b)->arr[2])

#define vector_equal(a,b) equal((a)->arr[0], (b)->arr[0]) && equal((a)->arr[1], (b)->arr[1]) && equal((a)->arr[2], (b)->arr[2])

#define matrix_equal(a,b) equal((a)->arr[0], (b)->arr[0])\
    && equal((a)->arr[1], (b)->arr[1])\
    && equal((a)->arr[2], (b)->arr[2])\
    && equal((a)->arr[3], (b)->arr[3])\
    && equal((a)->arr[4], (b)->arr[4])\
    && equal((a)->arr[5], (b)->arr[5])\
    && equal((a)->arr[6], (b)->arr[6])\
    && equal((a)->arr[7], (b)->arr[7])\
    && equal((a)->arr[8], (b)->arr[8])\
    && equal((a)->arr[9], (b)->arr[9])\
    && equal((a)->arr[10], (b)->arr[10])\
    && equal((a)->arr[11], (b)->arr[11])\
    && equal((a)->arr[12], (b)->arr[12])\
    && equal((a)->arr[13], (b)->arr[13])\
    && equal((a)->arr[14], (b)->arr[14])\
    && equal((a)->arr[15], (b)->arr[15])

typedef struct pt {
    double arr[4];
} *Point;

// TODO implement construction/destruction
typedef struct pts {
    Point points;
    size_t points_num;
} *Points;

typedef struct v {
    double arr[4];
} *Vector;

/*
 * 0  1  2  3
 * 4  5  6  7
 * 8  9  A  B
 * C  D  E  F
 *
 * 00 01 02 03
 * 10 11 12 13
 * 20 21 22 23
 * 30 31 32 33
 *
 * 0,0 -> 0 = 0 * 4 + 0
 * 0,1 -> 1 = 0 * 4 + 1
 * 0,2 -> 2 = 0 * 4 + 2
 * 0,3 -> 3 = 0 * 4 + 3
 * 1,0 -> 4 = 1 * 4 + 0
 * 1,1 -> 5 = 1 * 4 + 1
 * 1,2 -> 6 = 1 * 4 + 2
 * 1,3 -> 7 = 1 * 4 + 3
 * 2,0 -> 8 = 2 * 4 + 0
 * 2,1 -> 9 = 2 * 4 + 1
 * 2,2 -> A = 2 * 4 + 2
 * 2,3 -> B = 2 * 4 + 3
 * 3,0 -> C = 3 * 4 + 0
 * 3,1 -> D = 3 * 4 + 1
 * 3,2 -> E = 3 * 4 + 2
 * 3,3 -> F = 3 * 4 + 3
 */
typedef struct m {
    double arr[16];
} *Matrix;


// utilities
#define linalg_null_check_void(x)\
    if ((x) == NULL || (x)->arr == NULL) {\
        printf("Found a NULL pointer.\n");\
        return;\
    }

#define linalg_null_check(x,rv)\
    if ((x) == NULL || (x)->arr == NULL) {\
        printf("Found a NULL pointer.\n");\
        return (rv);\
    }

static const double MATRIX_IDENTITY[16] = {
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.0,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0
};

int point_to_string(char *s, size_t n, Point pt);

int vector_to_string(char *s, size_t n, Vector v);

int matrix_to_string(char *s, size_t n, Matrix m);


// construction
Point point(double x, double y, double z);
Point point_copy(Point pt);
Point point_default();

Vector vector(double x, double y, double z);
Vector vector_copy(Vector v);
Vector vector_default();

void matrix(double aa, double ab, double ac, double ad,
              double ba, double bb, double bc, double bd,
              double ca, double cb, double cc, double cd,
              double da, double db, double dc, double dd,
              Matrix res);

Matrix matrix_alloc(double aa, double ab, double ac, double ad,
              double ba, double bb, double bc, double bd,
              double ca, double cb, double cc, double cd,
              double da, double db, double dc, double dd);

void matrix_copy(Matrix m, Matrix res);
Matrix matrix_copy_alloc(Matrix m);
Matrix matrix_default();

// destruction
void point_free(Point pt);

void vector_free(Vector v);

void matrix_free(Matrix m);


// subtract pt2 from pt1
void vector_from_points(Point pt1, Point pt2, Vector res);
Vector vector_from_points_alloc(Point pt1, Point pt2);

double vector_magnitude(Vector v);

void vector_normalize(Vector v, Vector res);
Vector vector_normalize_alloc(Vector v);

double vector_dot(Vector a, Vector b);

void vector_cross(Vector a, Vector b, Vector res);
Vector vector_cross_alloc(Vector a, Vector b);

void vector_reflect(Vector a, Vector b, Vector res);
Vector vector_reflect_alloc(Vector a, Vector b);

void vector_scale(Vector input, double scalar);

void matrix_identity(Matrix res);
Matrix matrix_identity_alloc();

void matrix_translate(double x, double y, double z, Matrix res);
Matrix matrix_translate_alloc(double x, double y, double z);

void matrix_scale(double x, double y, double z, Matrix res);
Matrix matrix_scale_alloc(double x, double y, double z);

void matrix_rotate_x(double rad, Matrix res);
Matrix matrix_rotate_x_alloc(double rad);

void matrix_rotate_y(double rad, Matrix res);
Matrix matrix_rotate_y_alloc(double rad);

void matrix_rotate_z(double rad, Matrix res);
Matrix matrix_rotate_z_alloc(double rad);

void matrix_shear(double xy, double xz, double yx, double yz, double zx, double zy, Matrix res);
Matrix matrix_shear_alloc(double xy, double xz, double yx, double yz, double zx, double zy);

void matrix_multiply(Matrix a, Matrix b, Matrix res);
Matrix matrix_multiply_alloc(Matrix a, Matrix b);

void matrix_vector_multiply(Matrix a, Vector b, Vector res);
Vector matrix_vector_multiply_alloc(Matrix a, Vector b);

void matrix_point_multiply(Matrix a, Point b, Point res);
Point matrix_point_multiply_alloc(Matrix a, Point b);

void matrix_transpose(Matrix m, Matrix res);
Matrix matrix_transpose_alloc(Matrix m);

void matrix_inverse(Matrix m, Matrix res);
Matrix matrix_inverse_alloc(Matrix m);
#endif
