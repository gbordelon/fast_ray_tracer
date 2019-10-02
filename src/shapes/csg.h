#ifndef CSG
#define CSG

#include "shapes.h"

Shape csg_alloc(enum csg_ops_enum op, Shape left_child, Shape right_child);
void csg(Shape s, enum csg_ops_enum op, Shape left_child, Shape right_child);

#endif
