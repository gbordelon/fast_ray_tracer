#ifndef OBJ_LOADER
#define OBJ_LOADER

#include "../../shapes/shapes.h"

void construct_group_from_obj_file(const char * file_path, bool use_mtl, Shape result_group);

#endif
