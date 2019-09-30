#ifndef LINALG
#define LINALG

#include <stdio.h>
#include <stddef.h>

#define EPSILON 0.00001
#define equal(a,b) (fabs((a) - (b)) < EPSILON)

#define point_equal(a,b) equal((a)[0], (b)[0]) && equal((a)[1], (b)[1]) && equal((a)[2], (b)[2])

#define vector_equal(a,b) equal((a)[0], (b)[0]) && equal((a)[1], (b)[1]) && equal((a)[2], (b)[2])

#define matrix_equal(a,b) equal((a)[0], (b)[0])\
    && equal((a)[1], (b)[1])\
    && equal((a)[2], (b)[2])\
    && equal((a)[3], (b)[3])\
    && equal((a)[4], (b)[4])\
    && equal((a)[5], (b)[5])\
    && equal((a)[6], (b)[6])\
    && equal((a)[7], (b)[7])\
    && equal((a)[8], (b)[8])\
    && equal((a)[9], (b)[9])\
    && equal((a)[10], (b)[10])\
    && equal((a)[11], (b)[11])\
    && equal((a)[12], (b)[12])\
    && equal((a)[13], (b)[13])\
    && equal((a)[14], (b)[14])\
    && equal((a)[15], (b)[15])

typedef double Point[4];
typedef double Vector[4];
typedef double Matrix[16];


typedef struct pts {
    Point *points;
    size_t points_num;
} *Points;

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


static const double POINT_IDENTITY[4] = {
    0.0,0.0,0.0,1.0
};

static const double VECTOR_IDENTITY[4] = {
    0.0,0.0,0.0,0.0
};

static const double MATRIX_IDENTITY[16] = {
    1.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.0,
    0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,1.0
};

//#define matrix_identity() {1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0}
#define matrix_identity(m) memcpy((m), MATRIX_IDENTITY, sizeof(Matrix))
#define point_default(p) memcpy((p), POINT_IDENTITY, sizeof(Point))
#define vector_default(v) memcpy((v), VECTOR_IDENTITY, sizeof(Vector))

#define point(x,y,z,res) (res)[0]=(x);(res)[1]=(y);(res)[2]=(z);(res)[3]=1.0
#define vector(x,y,z,res) (res)[0]=(x);(res)[1]=(y);(res)[2]=(z);(res)[3]=0.0

void point_print(Point p);
void vector_print(Vector v);
void matrix_print(Matrix m);


// construction
void point_copy(Point to, Point from);

void vector_copy(Vector to, Vector from);

void matrix(double aa, double ab, double ac, double ad,
              double ba, double bb, double bc, double bd,
              double ca, double cb, double cc, double cd,
              double da, double db, double dc, double dd,
              Matrix res);

void matrix_copy(const Matrix m, Matrix res);


// subtract pt2 from pt1
void vector_from_points(Point pt1, Point pt2, Vector res);

void vector_cross(Vector a, Vector b, Vector res);

double vector_magnitude(Vector v);

void vector_normalize(Vector v, Vector res);

double vector_dot(Vector a, Vector b);

void vector_cross(Vector a, Vector b, Vector res);

void vector_reflect(Vector a, Vector b, Vector res);

void vector_scale(Vector input, double scalar);

void matrix_translate(double x, double y, double z, Matrix res);

void matrix_scale(double x, double y, double z, Matrix res);

void matrix_rotate_x(double rad, Matrix res);

void matrix_rotate_y(double rad, Matrix res);

void matrix_rotate_z(double rad, Matrix res);

void matrix_shear(double xy, double xz, double yx, double yz, double zx, double zy, Matrix res);

void matrix_multiply(const Matrix a, const Matrix b, Matrix res);
void transform_chain(const Matrix a, Matrix b);

void matrix_vector_multiply(const Matrix a, const Vector b, Vector res);

void matrix_point_multiply(const Matrix a, const Point b, Point res);

void matrix_transpose(const Matrix m, Matrix res);

void matrix_inverse(const Matrix m, Matrix res);

#endif
