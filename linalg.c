#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "linalg.h"



void
point_copy(Point to, Point from)
{
    memcpy(to, from, sizeof(Point));
}

void
vector_copy(Vector to, Vector from)
{
    memcpy(to, from, sizeof(Vector));
}

void
vector_from_points(Point pt1, Point pt2, Vector res)
{
    res[0] = pt1[0] - pt2[0];
    res[1] = pt1[1] - pt2[1];
    res[2] = pt1[2] - pt2[2];
    res[3] = 0.0;
}

void
vector_reflect(Vector input, Vector normal, Vector res)
{
    double ddot = 2 * vector_dot(input, normal);
    res[0] = input[0] - normal[0] * ddot;
    res[1] = input[1] - normal[1] * ddot;
    res[2] = input[2] - normal[2] * ddot;
    res[3] = 0.0;
}

void
vector_scale(Vector input, double scalar)
{
    input[0] *= scalar;
    input[1] *= scalar;
    input[2] *= scalar;
}

void
matrix(double aa, double ab, double ac, double ad,
       double ba, double bb, double bc, double bd,
       double ca, double cb, double cc, double cd,
       double da, double db, double dc, double dd,
       Matrix res)
{
    res[0] = aa;
    res[1] = ab;
    res[2] = ac;
    res[3] = ad;

    res[4] = ba;
    res[5] = bb;
    res[6] = bc;
    res[7] = bd;

    res[8] = ca;
    res[9] = cb;
    res[10] = cc;
    res[11] = cd;

    res[12] = da;
    res[13] = db;
    res[14] = dc;
    res[15] = dd;
}

void
matrix_copy(const Matrix m, Matrix res)
{
    memcpy(res, m, sizeof(Matrix));
}

int
array_to_string(char *buf, size_t n, const char *name, double arr[4])
{
    return snprintf(buf, n, "%s: [%f %f %f %f]", name, arr[0], arr[1], arr[2], arr[3]);
}

int
matrix_to_string(char *buf, size_t n, Matrix m)
{
    return snprintf(buf, n, "Matrix: [\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n %f %f %f %f\n]",
                    m[0], m[1], m[2], m[3],
                    m[4], m[5], m[6], m[7],
                    m[8], m[9], m[10], m[11],
                    m[12], m[13], m[14], m[15]);
}

void
array_print(double arr[4], const char *name)
{
    static char buf[256];
    array_to_string(buf, 256, name, arr);
    printf("%s\n",buf);
}

void
vector_print(Vector v)
{
    array_print(v, "Vector");
}

void
point_print(Vector v)
{
    array_print(v, "Point");
}

void
matrix_print(Matrix m)
{
    static char buf[256];
    matrix_to_string(buf, 256, m);
    printf("%s\n",buf);
}

double
vector_magnitude(Vector v)
{
    double sum;

    if (v[3] != 0.0) {
        // error
    }

    sum = v[0] * v[0] +
          v[1] * v[1] +
          v[2] * v[2];
    return sqrt(sum);
}

void
vector_normalize(Vector v, Vector res)
{
    double inv = 1.0 / vector_magnitude(v);
    vector_copy(res, v);
    vector_scale(res, inv);
    res[3] = 0.0;
}

double
vector_dot(Vector a, Vector b)
{
    return a[0] * b[0] +
           a[1] * b[1] +
           a[2] * b[2];
}

void
vector_cross(Vector a, Vector b, Vector res)
{
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
    res[3] = 0.0;
}

void
matrix_translate(double x, double y, double z, Matrix res)
{
    matrix_identity(res);
    res[3] = x;
    res[7] = y;
    res[11] = z;
}

void
matrix_scale(double x, double y, double z, Matrix res)
{
    matrix_identity(res);
    res[0] = x;
    res[5] = y;
    res[10] = z;
}

void
matrix_rotate_x(double rad, Matrix res)
{
    matrix_identity(res);
    double cr = cos(rad);
    res[5] = cr;
    res[10] = cr;
    res[6] = -sin(rad);
    res[9] = sin(rad);
}

void
matrix_rotate_y(double rad, Matrix res)
{
    matrix_identity(res);
    double cr = cos(rad);
    res[0] = cr;
    res[10] = cr;
    res[8] = -sin(rad);
    res[2] = sin(rad);
}

void
matrix_rotate_z(double rad, Matrix res)
{
    matrix_identity(res);
    double cr = cos(rad);
    res[0] = cr;
    res[5] = cr;
    res[1] = -sin(rad);
    res[4] = sin(rad);
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
    res[1] = xy;
    res[2] = xz;
    res[4] = yx;
    res[6] = yz;
    res[8] = zx;
    res[9] = zy;
}

void
matrix_multiply(const Matrix a, const Matrix b, Matrix res)
{
    int i, j;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            res[i*4+j] = a[i*4+0] * b[0*4+j] +
                              a[i*4+1] * b[1*4+j] +
                              a[i*4+2] * b[2*4+j] +
                              a[i*4+3] * b[3*4+j];
        }
    }
}

void
transform_chain(const Matrix a, Matrix b)
{
    Matrix tmp;
    matrix_multiply(a,b,tmp);
    matrix_copy(tmp,b);
}

void
matrix_array_multiply(const Matrix a, const double b[4], double res[4])
{
    int i;

    for (i = 0; i < 4; i++) {
        res[i] = a[i*4+0] * b[0] +
                 a[i*4+1] * b[1] +
                 a[i*4+2] * b[2] +
                 a[i*4+3] * b[3];
    }
}

void
matrix_point_multiply(const Matrix a, const Point b, Point res)
{
    matrix_array_multiply(a, b, res);
}

void
matrix_vector_multiply(const Matrix a, const Vector b, Vector res)
{
    matrix_array_multiply(a, b, res);
}

void
matrix_transpose(const Matrix m, Matrix res)
{
    res[0]  = m[0];
    res[1]  = m[4];
    res[2]  = m[8];
    res[3]  = m[12];
    res[4]  = m[1];
    res[5]  = m[5];
    res[6]  = m[9];
    res[7]  = m[13];
    res[8]  = m[2];
    res[9]  = m[6];
    res[10] = m[10];
    res[11] = m[14];
    res[12] = m[3];
    res[13] = m[7];
    res[14] = m[11];
    res[15] = m[15];
}

void
matrix_inverse(const Matrix m, Matrix res)
{
    double det;

    double sub05 = m[10] * m[15] - m[11] * m[14];
    double sub06 = m[9] * m[15] - m[11] * m[13];
    double sub07 = m[9] * m[14] - m[10] * m[13];
    double sub16 = m[8] * m[15] - m[11] * m[12];
    double sub17 = m[8] * m[14] - m[10] * m[12];
    double sub27 = m[8] * m[13] - m[9] * m[12];

    // minor of 00
    double sub0 = m[5] * sub05 - m[6] * sub06 + m[7] * sub07;
    // minor of 01
    double sub1 = m[4] * sub05 - m[6] * sub16 + m[7] * sub17;
    // minor of 02
    double sub2 = m[4] * sub06 - m[5] * sub16 + m[7] * sub27;
    // minor of 03
    double sub3 = m[4] * sub07 - m[5] * sub17 + m[6] * sub27;


    // minor 10
    double sub4 = m[1] * sub05 - m[2] * sub06 + m[3] * sub07;
    // minor 11
    double sub5 = m[0] * sub05 - m[2] * sub16 + m[3] * sub17;
    // minor 12
    double sub6 = m[0] * sub06 - m[1] * sub16 + m[3] * sub27;
    // minor 13
    double sub7 = m[0] * sub07 - m[1] * sub17 + m[2] * sub27;


    double sub81 = m[6] * m[15] - m[7] * m[14];
    double sub82 = m[5] * m[15] - m[7] * m[13];
    double sub83 = m[5] * m[14] - m[6] * m[13];
    double sub92 = m[4] * m[15] - m[7] * m[12];
    double sub93 = m[4] * m[14] - m[6] * m[12];
    double subA3 = m[4] * m[13] - m[5] * m[12];

    // minor 20
    double sub8 = m[1] * sub81 - m[2] * sub82 + m[3] * sub83;
    // minor 21
    double sub9 = m[0] * sub81 - m[2] * sub92 + m[3] * sub93;
    // minor 22
    double subA = m[0] * sub82 - m[1] * sub92 + m[3] * subA3;
    // minor 23
    double subB = m[0] * sub83 - m[1] * sub93 + m[2] * subA3;


    double subC1 = m[6] * m[11] - m[7] * m[10];
    double subC2 = m[5] * m[11] - m[7] * m[9];
    double subC3 = m[5] * m[10] - m[6] * m[9];
    double subD2 = m[4] * m[11] - m[7] * m[8];
    double subD3 = m[4] * m[10] - m[6] * m[8];
    double subE3 = m[4] * m[9] - m[5] * m[8];

    // minor30
    double subC = m[1] * subC1 - m[2] * subC2 + m[3] * subC3;
    // minor31
    double subD = m[0] * subC1 - m[2] * subD2 + m[3] * subD3;
    // minor32
    double subE = m[0] * subC2 - m[1] * subD2 + m[3] * subE3;
    // minor33
    double subF = m[0] * subC3 - m[1] * subD3 + m[2] * subE3;

    det =  m[0] * sub0 - m[1] * sub1 + m[2] * sub2 - m[3] * sub3;

    if (equal(det,0.0)) {
        printf("determinant is zero\n");
    }

    res[0] = sub0 / det;
    res[4] = -sub1 / det;
    res[8] = sub2 / det;
    res[12] = -sub3 / det;
    res[1] = -sub4 / det;
    res[5] = sub5 / det;
    res[9] = -sub6 / det;
    res[13] = sub7 / det;
    res[2] = sub8 / det;
    res[6] = -sub9 / det;
    res[10] = subA / det;
    res[14] = -subB / det;
    res[3] = -subC / det;
    res[7] = subD / det;
    res[11] = -subE / det;
    res[15] = subF / det;
}
