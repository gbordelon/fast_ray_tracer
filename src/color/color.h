#ifndef COLOR
#define COLOR

typedef double Color[4];

static const Color BLACK = {
    0.0,
    0.0,
    0.0,
    0.0
};

static const Color WHITE = {
    1.0,
    1.0,
    1.0,
    0.0
};

#define color_default(c) memcpy((c), BLACK, sizeof(Color));
#define color(r,g,b) { (r), (g), (b), 0.0 }

void color_accumulate(Color acc, const Color other);
void color_scale(Color acc, const double scalar);
void color_copy(Color to, const Color from);

void print_color(const Color c);

void rgb_to_lab(const Color rgb, Color lab);
void lab_to_rgb(const Color lab, Color rgb);

int lab_compare_l(const Color l, const Color r);
int lab_compare_a(const Color l, const Color r);
int lab_compare_b(const Color l, const Color r);

#endif
