#ifndef RAY
#define RAY

typedef struct ray {
    Point origin;
    Vector direction;
} *Ray;

void ray_array(Point origin, Vector direction, Ray ray);
int ray_to_string(char *s, size_t n, Ray r);
void ray_transform(Ray original, Matrix m, Ray res);
void ray_position(Ray ray, double t, Point position);

#endif
