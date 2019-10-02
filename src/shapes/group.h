#ifndef GROUP
#define GROUP

#include "shapes.h"

Shape group_alloc(Shape children, size_t num_children);
void group(Shape s, Shape children, size_t num_children);
void group_add_children(Shape group, Shape children, size_t num_children);

#endif
