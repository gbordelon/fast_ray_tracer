#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "linalg.h"



/*
 * Create a data structure which will store n Matrixes when the first matrix is alloc'd
 * When it gets full, double the size and memcpy over to the new block
 *      all matrix pointers will be broken at this point
 *      instead of using matrix pointers, use an index?
 *          this fixes the broken pointer problem but how do i defrag?
 * keep a free list of matrixes which have been released
 * use free list first when a new matrix is being alloc'd
 * 
 * Will this actually be faster than malloc and free?
 */


Point
point_default()
{
    Point pt = (Point) malloc(sizeof(struct pt));
    if (pt == NULL) {
        // error check
    }

    pt->arr[3] = 1.0;
    return pt;
}

Point
point(double x,
      double y,
      double z)
{
    Point pt = point_default();
    linalg_null_check(pt,NULL)

    pt->arr[0] = x;
    pt->arr[1] = y;
    pt->arr[2] = z;

    return pt;
}

Point
point_copy(Point pt)
{
    linalg_null_check(pt,NULL)
    return point(pt->arr[0], pt->arr[1], pt->arr[2]);
}

Vector
vector_default()
{
    Vector v = (Vector) malloc(sizeof(struct v));
    if (v == NULL) {
        // error check
    }

    v->arr[3] = 0.0;
    return v;
}

Vector
vector(double x,
       double y,
       double z)
{
    Vector v = vector_default();
    linalg_null_check(v,NULL)

    v->arr[0] = x;
    v->arr[1] = y;
    v->arr[2] = z;

    return v;
}

Vector
vector_copy(Vector v)
{
    linalg_null_check(v,NULL)
    return vector(v->arr[0], v->arr[1], v->arr[2]);
}

void
vector_from_points(Point pt1, Point pt2, Vector res)
{
    linalg_null_check_void(pt1)
    linalg_null_check_void(pt2)
    linalg_null_check_void(res)

    res->arr[0] = pt1->arr[0] - pt2->arr[0];
    res->arr[1] = pt1->arr[1] - pt2->arr[1];
    res->arr[2] = pt1->arr[2] - pt2->arr[2];
}

Vector
vector_from_points_alloc(Point pt1, Point pt2)
{
    linalg_null_check(pt1,NULL)
    linalg_null_check(pt2,NULL)

    return vector_from_arrays_alloc(pt1->arr, pt2->arr);
}

Vector
vector_from_arrays_alloc(double pt1[4], double pt2[4])
{
    Vector v = vector_default();

    linalg_null_check(v,NULL)

    v->arr[0] = pt1[0] - pt2[0];
    v->arr[1] = pt1[1] - pt2[1];
    v->arr[2] = pt1[2] - pt2[2];
    // v->arr[3] is already set

    return v;
}

void
array_from_arrays(double pt1[4], double pt2[4], double res[4])
{
    res[0] = pt1[0] - pt2[0];
    res[1] = pt1[1] - pt2[1];
    res[2] = pt1[2] - pt2[2];
}

void
vector_reflect(Vector input, Vector normal, Vector res)
{
    double ddot = 2 * vector_dot(input, normal);
    res->arr[0] = input->arr[0] - normal->arr[0] * ddot;
    res->arr[1] = input->arr[1] - normal->arr[1] * ddot;
    res->arr[2] = input->arr[2] - normal->arr[2] * ddot;
}

Vector
vector_reflect_alloc(Vector a, Vector b)
{
    Vector res = vector_default();
    vector_reflect(a, b, res);
    return res;
}

void
vector_scale(Vector input, double scalar)
{
    input->arr[0] *= scalar;
    input->arr[1] *= scalar;
    input->arr[2] *= scalar;
}

Matrix
matrix_default()
{
    Matrix m = (Matrix) malloc(sizeof(struct m));
    if (m == NULL) {
        // error check
    }

    return m;
}

void
matrix(double aa, double ab, double ac, double ad,
       double ba, double bb, double bc, double bd,
       double ca, double cb, double cc, double cd,
       double da, double db, double dc, double dd,
       Matrix res)
{
    linalg_null_check_void(res)

    res->arr[0] = aa;
    res->arr[1] = ab;
    res->arr[2] = ac;
    res->arr[3] = ad;

    res->arr[4] = ba;
    res->arr[5] = bb;
    res->arr[6] = bc;
    res->arr[7] = bd;

    res->arr[8] = ca;
    res->arr[9] = cb;
    res->arr[10] = cc;
    res->arr[11] = cd;

    res->arr[12] = da;
    res->arr[13] = db;
    res->arr[14] = dc;
    res->arr[15] = dd;
}

Matrix
matrix_alloc(double aa, double ab, double ac, double ad,
       double ba, double bb, double bc, double bd,
       double ca, double cb, double cc, double cd,
       double da, double db, double dc, double dd)
{
    Matrix res = matrix_default();
    matrix(aa, ab, ac, ad,
           ba, bb, bc, bd,
           ca, cb, cc, cd,
           da, db, dc, dd,
           res);
    return res;
}

void
matrix_copy(Matrix m, Matrix res)
{
    linalg_null_check_void(m)
    linalg_null_check_void(res)

    memcpy(res->arr, m->arr, 16 * sizeof(double));
}

Matrix
matrix_copy_alloc(Matrix m)
{
    Matrix res = matrix_default();
    matrix_copy(m, res);
    return res;
}

void
point_free(Point pt)
{
    if (pt != NULL) {
        free(pt);
    }
}

void
vector_free(Vector v)
{
    if (v != NULL) {
        free(v);
    }
}

void
matrix_free(Matrix m)
{
    if (m != NULL) {
        free(m);
    }
}

int
point_to_string(char * buf, size_t n, Point pt)
{
    linalg_null_check(pt,0)
    return snprintf(buf, n, "Point: [%f %f %f]", pt->arr[0], pt->arr[1], pt->arr[2]);
}

int
vector_to_string(char * buf, size_t n, Vector v)
{
    linalg_null_check(v,0)
    return snprintf(buf, n, "Vector: [%f %f %f]", v->arr[0], v->arr[1], v->arr[2]);
}

int
matrix_to_string(char * buf, size_t n, Matrix m)
{
    linalg_null_check(m,0)
    return snprintf(buf, n, "Matrix: [\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n]",
                    m->arr[0], m->arr[1], m->arr[2], m->arr[3],
                    m->arr[4], m->arr[5], m->arr[6], m->arr[7],
                    m->arr[8], m->arr[9], m->arr[10], m->arr[11],
                    m->arr[12], m->arr[13], m->arr[14], m->arr[15]);
}

double
vector_magnitude(Vector v)
{
    double sum;

    linalg_null_check(v, 0.0)

    if (v->arr[3] != 0.0) {
        // error
    }

    sum = v->arr[0] * v->arr[0] +
          v->arr[1] * v->arr[1] +
          v->arr[2] * v->arr[2];
    return sqrt(sum);
}

void
vector_normalize(Vector v, Vector res)
{
    double mag;
    int i;

    linalg_null_check_void(v)
    linalg_null_check_void(res)

    mag = vector_magnitude(v);
    for (i = 0; i < 3; i++) {
        res->arr[i] = v->arr[i] / mag;
    }
}

Vector
vector_normalize_alloc(Vector v)
{
    Vector retval = vector_default();
    vector_normalize(v, retval);
    return retval;
}

double
vector_dot(Vector a, Vector b)
{
    linalg_null_check(a,0.0)
    linalg_null_check(b,0.0)

    return array_dot(a->arr, b->arr);
}

double
array_dot(double a[4], double b[4])
{
    return a[0] * b[0] +
           a[1] * b[1] +
           a[2] * b[2];
}

void
vector_cross_arrays(double a[4], double b[4], Vector res)
{
    linalg_null_check_void(res)

    res->arr[0] = a[1] * b[2] - a[2] * b[1];
    res->arr[1] = a[2] * b[0] - a[0] * b[2];
    res->arr[2] = a[0] * b[1] - a[1] * b[0];
}

void
vector_cross(Vector a, Vector b, Vector res)
{
    vector_cross_arrays(a->arr, b->arr, res);
}

Vector
vector_cross_arrays_alloc(double a[4], double b[4])
{
    Vector v = vector_default();
    vector_cross_arrays(a, b, v);
    return v;
}

Vector
vector_cross_alloc(Vector a, Vector b)
{
    return vector_cross_arrays_alloc(a->arr,b->arr);
}

void
matrix_identity(Matrix res)
{
    linalg_null_check_void(res)
    memcpy(res->arr, MATRIX_IDENTITY, 16 * sizeof(double));
}

Matrix
matrix_identity_alloc()
{
    Matrix res = matrix_alloc(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
    linalg_null_check(res,NULL)
    return res;
}

void
matrix_translate(double x, double y, double z, Matrix res)
{
    matrix_identity(res);
    res->arr[3] = x;
    res->arr[7] = y;
    res->arr[11] = z;
}

Matrix
matrix_translate_alloc(double x, double y, double z)
{
    Matrix res = matrix_alloc(1,0,0,x, 0,1,0,y, 0,0,1,z, 0,0,0,1);
    linalg_null_check(res,NULL)
    return res;
}

void
matrix_scale(double x, double y, double z, Matrix res)
{
    matrix_identity(res);
    res->arr[0] = x;
    res->arr[5] = y;
    res->arr[10] = z;
}

Matrix
matrix_scale_alloc(double x, double y, double z)
{
    Matrix res = matrix_alloc(x,0,0,0, 0,y,0,0, 0,0,z,0, 0,0,0,1);
    linalg_null_check(res,NULL)
    return res;
}

void
matrix_rotate_x(double rad, Matrix res)
{
    matrix_identity(res);
    double cr = cos(rad);
    res->arr[5] = cr;
    res->arr[10] = cr;
    res->arr[6] = -sin(rad);
    res->arr[9] = sin(rad);
}

Matrix
matrix_rotate_x_alloc(double rad)
{
    double cr = cos(rad);
    Matrix res = matrix_alloc(1,0,0,0,
                        0,cr,-sin(rad),0,
                        0,sin(rad),cr,0,
                        0,0,0,1);
    linalg_null_check(res,NULL)
    return res;
}

void
matrix_rotate_y(double rad, Matrix res)
{
    matrix_identity(res);
    double cr = cos(rad);
    res->arr[0] = cr;
    res->arr[10] = cr;
    res->arr[8] = -sin(rad);
    res->arr[2] = sin(rad);
}

Matrix
matrix_rotate_y_alloc(double rad)
{
    double cr = cos(rad);
    Matrix res = matrix_alloc(cr,0,sin(rad),0,
                        0,1,0,0,
                        -sin(rad),0,cr,0,
                        0,0,0,1);
    linalg_null_check(res,NULL)
    return res;
}

void
matrix_rotate_z(double rad, Matrix res)
{
    matrix_identity(res);
    double cr = cos(rad);
    res->arr[0] = cr;
    res->arr[5] = cr;
    res->arr[1] = -sin(rad);
    res->arr[4] = sin(rad);
}

Matrix
matrix_rotate_z_alloc(double rad)
{
    double cr = cos(rad);
    Matrix res = matrix_alloc(cr,-sin(rad),0,0,
                        sin(rad),cr,0,0,
                        0,0,1,0,
                        0,0,0,1);
    linalg_null_check(res,NULL)
    return res;
}

void
matrix_shear(double xy,
                  double xz,
                  double yx,
                  double yz,
                  double zx,
                  double zy,
                  Matrix res)
{
    matrix_identity(res);
    res->arr[1] = xy;
    res->arr[2] = xz;
    res->arr[4] = yx;
    res->arr[6] = yz;
    res->arr[8] = zx;
    res->arr[9] = zy;
}

Matrix
matrix_shear_alloc(double xy,
             double xz,
             double yx,
             double yz,
             double zx,
             double zy)
{
    Matrix res = matrix_alloc(1,xy,xz,0, yx,1,yz,0, zx,zy,1,0, 0,0,0,1);
    linalg_null_check(res,NULL)
    return res;
}

void
matrix_multiply(Matrix a, Matrix b, Matrix res)
{
    int i, j;

    linalg_null_check_void(a)
    linalg_null_check_void(b)
    linalg_null_check_void(res)

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            res->arr[i*4+j] = a->arr[i*4+0] * b->arr[0*4+j] +
                              a->arr[i*4+1] * b->arr[1*4+j] +
                              a->arr[i*4+2] * b->arr[2*4+j] +
                              a->arr[i*4+3] * b->arr[3*4+j];
        }
    }
}

Matrix
matrix_multiply_alloc(Matrix a, Matrix b)
{
    Matrix res = matrix_default();
    matrix_multiply(a,b,res);
    return res;
}

void
matrix_point_multiply(Matrix a, Point b, Point res)
{
    linalg_null_check_void(a)
    linalg_null_check_void(b)
    linalg_null_check_void(res)

    matrix_array_multiply(a, b->arr, res->arr);
}

Point
matrix_point_multiply_alloc(Matrix a, Point b)
{
    Point res = point_default();
    matrix_point_multiply(a,b,res);
    return res;
}

void
matrix_vector_multiply(Matrix a, Vector b, Vector res)
{
    linalg_null_check_void(a)
    linalg_null_check_void(b)
    linalg_null_check_void(res)

    matrix_array_multiply(a, b->arr, res->arr);
}

Vector
matrix_vector_multiply_alloc(Matrix a, Vector b)
{
    Vector res = vector_default();
    matrix_vector_multiply(a,b,res);
    return res;
}

void
matrix_array_multiply(Matrix a, double b[4], double res[4])
{
    int i;

    for (i = 0; i < 4; i++) {
        res[i] = a->arr[i*4+0] * b[0] +
                 a->arr[i*4+1] * b[1] +
                 a->arr[i*4+2] * b[2] +
                 a->arr[i*4+3] * b[3];
    }
}

void
matrix_transpose(Matrix m, Matrix res)
{
    linalg_null_check_void(m)
    linalg_null_check_void(res)

    res->arr[0]  = m->arr[0];
    res->arr[1]  = m->arr[4];
    res->arr[2]  = m->arr[8];
    res->arr[3]  = m->arr[12];
    res->arr[4]  = m->arr[1];
    res->arr[5]  = m->arr[5];
    res->arr[6]  = m->arr[9];
    res->arr[7]  = m->arr[13];
    res->arr[8]  = m->arr[2];
    res->arr[9]  = m->arr[6];
    res->arr[10] = m->arr[10];
    res->arr[11] = m->arr[14];
    res->arr[12] = m->arr[3];
    res->arr[13] = m->arr[7];
    res->arr[14] = m->arr[11];
    res->arr[15] = m->arr[15];
}

Matrix
matrix_transpose_alloc(Matrix m)
{
    Matrix res = matrix_default();
    matrix_transpose(m, res);
    return res;
}

void
matrix_inverse(Matrix m, Matrix res)
{
    double det;

    linalg_null_check_void(m)
    linalg_null_check_void(res)

    double *arr = m->arr;

    double sub05 = arr[10] * arr[15] - arr[11] * arr[14];
    double sub06 = arr[9] * arr[15] - arr[11] * arr[13];
    double sub07 = arr[9] * arr[14] - arr[10] * arr[13];
    double sub16 = arr[8] * arr[15] - arr[11] * arr[12];
    double sub17 = arr[8] * arr[14] - arr[10] * arr[12];
    double sub27 = arr[8] * arr[13] - arr[9] * arr[12];

    // minor of 00
    double sub0 = arr[5] * sub05 - arr[6] * sub06 + arr[7] * sub07;
    // minor of 01
    double sub1 = arr[4] * sub05 - arr[6] * sub16 + arr[7] * sub17;
    // minor of 02
    double sub2 = arr[4] * sub06 - arr[5] * sub16 + arr[7] * sub27;
    // minor of 03
    double sub3 = arr[4] * sub07 - arr[5] * sub17 + arr[6] * sub27;


    // minor 10
    double sub4 = arr[1] * sub05 - arr[2] * sub06 + arr[3] * sub07;
    // minor 11
    double sub5 = arr[0] * sub05 - arr[2] * sub16 + arr[3] * sub17;
    // minor 12
    double sub6 = arr[0] * sub06 - arr[1] * sub16 + arr[3] * sub27;
    // minor 13
    double sub7 = arr[0] * sub07 - arr[1] * sub17 + arr[2] * sub27;


    double sub81 = arr[6] * arr[15] - arr[7] * arr[14];
    double sub82 = arr[5] * arr[15] - arr[7] * arr[13];
    double sub83 = arr[5] * arr[14] - arr[6] * arr[13];
    double sub92 = arr[4] * arr[15] - arr[7] * arr[12];
    double sub93 = arr[4] * arr[14] - arr[6] * arr[12];
    double subA3 = arr[4] * arr[13] - arr[5] * arr[12];

    // minor 20
    double sub8 = arr[1] * sub81 - arr[2] * sub82 + arr[3] * sub83;
    // minor 21
    double sub9 = arr[0] * sub81 - arr[2] * sub92 + arr[3] * sub93;
    // minor 22
    double subA = arr[0] * sub82 - arr[1] * sub92 + arr[3] * subA3;
    // minor 23
    double subB = arr[0] * sub83 - arr[1] * sub93 + arr[2] * subA3;


    double subC1 = arr[6] * arr[11] - arr[7] * arr[10];
    double subC2 = arr[5] * arr[11] - arr[7] * arr[9];
    double subC3 = arr[5] * arr[10] - arr[6] * arr[9];
    double subD2 = arr[4] * arr[11] - arr[7] * arr[8];
    double subD3 = arr[4] * arr[10] - arr[6] * arr[8];
    double subE3 = arr[4] * arr[9] - arr[5] * arr[8];

    // minor30
    double subC = arr[1] * subC1 - arr[2] * subC2 + arr[3] * subC3;
    // minor31
    double subD = arr[0] * subC1 - arr[2] * subD2 + arr[3] * subD3;
    // minor32
    double subE = arr[0] * subC2 - arr[1] * subD2 + arr[3] * subE3;
    // minor33
    double subF = arr[0] * subC3 - arr[1] * subD3 + arr[2] * subE3;

    det =  arr[0] * sub0 - arr[1] * sub1 + arr[2] * sub2 - arr[3] * sub3;

    if (equal(det,0.0)) {
        printf("determinant is zero\n");
        // error
    }

    res->arr[0] = sub0 / det;
    res->arr[4] = -sub1 / det;
    res->arr[8] = sub2 / det;
    res->arr[12] = -sub3 / det;
    res->arr[1] = -sub4 / det;
    res->arr[5] = sub5 / det;
    res->arr[9] = -sub6 / det;
    res->arr[13] = sub7 / det;
    res->arr[2] = sub8 / det;
    res->arr[6] = -sub9 / det;
    res->arr[10] = subA / det;
    res->arr[14] = -subB / det;
    res->arr[3] = -subC / det;
    res->arr[7] = subD / det;
    res->arr[11] = -subE / det;
    res->arr[15] = subF / det;
}

Matrix
matrix_inverse_alloc(Matrix m)
{
    Matrix res = matrix_default();
    matrix_inverse(m, res);
    return res;
}
