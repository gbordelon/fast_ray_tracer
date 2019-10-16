#ifndef OBJ_LOADER
#define OBJ_LOADER

#include "../../color/color.h"
#include "../../shapes/shapes.h"

void construct_group_from_obj_file(const char * file_path, void (*color_space_fn)(const Color, Color), Shape result_group);

#endif
