from renderer import Camera, Light, allocate_lights
from shapes import Shape, Group, Material, allocate_shapes
from config import GlobalConfig
from transform import Transform

from copy import deepcopy
import yaml
from yaml import CLoader as Loader

def yaml_file_to_world_objects(file_path):
    tree = None
    with open(file_path, 'r') as f:
        tree = yaml.load(f, Loader=Loader)

    if tree is None:
        return []

    rv = {'camera':None,
          'lights':[],
          'world':[],
          'config':None}

    defines = {}
    extends_map = {}

    for obj in tree:
        if "define" in obj:
            k = obj["define"]
            v = obj.get("value")
            opt = obj.get("extend")
            defines[k] = v
            if opt is not None:
                extends_map[k] = opt

    # replace 'extends' in defines map
    for obj_name in extends_map:
        parent_name = extends_map[obj_name] # name of object which will be extended
        parent_value = defines[parent_name]
        child_value = defines[obj_name] # name of object with 'extends' keyword
        new_parent_value = deepcopy(parent_value)
        if type(new_parent_value) == dict:
            # assume child value is same type
            for k in child_value:
                new_parent_value[k] = child_value[k]
            defines[obj_name] = new_parent_value

    expand_defines_in_tree(tree, defines)

    for obj in tree:
        if "add" in obj:
            if obj["add"] == "camera":
                rv['camera'] = Camera.from_yaml(obj)
            elif obj["add"] == "light":
                rv['lights'].append(Light.from_yaml(obj))
            elif obj['add'] == 'config':
                rv['config'] = GlobalConfig.from_yaml(obj)
            else:
                possible_item = recursive_add(obj, defines)
                if possible_item is not None:
                    rv['world'].append(possible_item)

    return rv

def recursive_add(tree, defines):
    return Shape.from_yaml(tree, defines)

# replace occurrences of previous defines in the tree
def expand_defines_in_tree(tree, defines):
    for obj in tree:
        for k in defines:
            if "value" in obj and k in obj["value"]:
                new_parent_value = deepcopy(defines[k])
                i = obj["value"].index(k)
                del(obj["value"][i])
                for item in new_parent_value:
                    obj["value"].insert(i, item)
                    i += 1
            if "material" in obj and k in obj["material"]:
                if type(obj["material"]) == str:
                    obj["material"] = deepcopy(defines[k])
                elif type(obj["material"]) == dict:
                    tmp = obj["material"]
                    obj["material"] = deepcopy(defines[k])
                    for j in tmp:
                        obj["material"][j] = tmp[j]
            if "transform" in obj and k in obj["transform"]:
                i = obj["transform"].index(k)
                del(obj["transform"][i])
                for item in deepcopy(defines[k]):
                    obj["transform"].insert(i, item)
                    i += 1
            if "extend" in obj and k in obj["extend"]:
                del(obj["extend"])
                if "value" in obj and k in obj["value"]:
                    obj["value"][k] = deepcopy(defines[k])
                elif "value" in obj:
                    if type(obj["value"]) == dict:
                        tmp = obj["value"]
                        obj["value"] = deepcopy(defines[k])
                        for j in tmp:
                            obj["value"][j] = tmp[j]
                    else:
                        for item in deepcopy(defines[k]):
                            obj["value"].insert(0, item)
                else:
                    obj["value"] = deepcopy(defines[k])
            if "add" in obj:
                if k == obj["add"]:
                    new_defines = deepcopy(defines[k])
                    if "add" in new_defines and new_defines["add"] == "group" and "children" in new_defines:
                        expand_defines_in_tree(new_defines["children"], defines)

                    if "add" in new_defines and new_defines["add"] == "csg" and "left" in new_defines:
                        expand_defines_in_tree([new_defines["left"]], defines)
                    if "add" in new_defines and new_defines["add"] == "csg" and "right" in new_defines:
                        expand_defines_in_tree([new_defines["right"]], defines)

                    for l in new_defines:
                        if l != "material" and l != "transform":
                            obj[l] = new_defines[l]
                        if l == "material" and "material" not in obj:
                            obj[l] = new_defines[l]
                        if l == "transform":
                            if "transform" not in obj:
                                obj[l] = new_defines[l]
                            else:
                                i = 0
                                for xform in new_defines[l]:
                                    obj[l].insert(i, xform)
                                    i += 1
                if obj["add"] == "group" and "children" in obj:
                    expand_defines_in_tree(obj["children"], defines)
                #if obj["add"] == "csg" and "left" in obj and "right" in obj:
                #    expand_defines_in_tree(obj["left"], defines)
                #    expand_defines_in_tree(obj["right"], defines)


def world_objects_to_c_file(obj):
    if 'config' not in obj or obj['config'] is None:
        obj['config'] = GlobalConfig.from_yaml(dict())
    
    c_code = """#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>

#include "src/libs/canvas/canvas.h"
#include "src/libs/linalg/linalg.h"
#include "src/libs/obj_loader/obj_loader.h"
#include "src/renderer/renderer.h"
#include "src/renderer/camera.h"
#include "src/renderer/world.h"
#include "src/shapes/shapes.h"
#include "src/shapes/cone.h"
#include "src/shapes/csg.h"
#include "src/shapes/cube.h"
#include "src/shapes/cylinder.h"
#include "src/shapes/group.h"
#include "src/shapes/plane.h"
#include "src/shapes/sphere.h"
#include "src/shapes/triangle.h"
#include "src/shapes/toroid.h"

int main()
{{
{0}
{1}
{2}
{4}
    Shape world_group = array_of_shapes(1);
    group(world_group, all_shapes, {5});
    world_group->divide(world_group, divide_threshold);

    World w = world();
    w->lights = all_lights;
    w->lights_num = {3};
    w->shapes = world_group;
    w->shapes_num = 1;
    Canvas c = render_multi(cam, w, cam->sample_num/*usteps*/, cam->sample_num/*vsteps*/, cam->aperture.jitter, thread_count);
    Ppm ppm = construct_ppm(c, true);

    FILE * pFile;
    pFile = fopen ("/tmp/{6}.ppm", "wb");
    fwrite (ppm->arr, sizeof(unsigned char), ppm->len, pFile);
    fclose (pFile);

    ppm_free(ppm);
    ppm = construct_ppm(c, false);

    pFile = fopen("/tmp/{7}.ppm", "wb");
    fwrite (ppm->arr, sizeof(unsigned char), ppm->len, pFile);
    fclose (pFile);

    ppm_free(ppm);
    canvas_free(c);

    return 0;
}}
""".format(obj['config'].c_repr(),
           obj['camera'].c_repr(),
           allocate_lights(obj['lights']),
           len(obj['lights']),
           allocate_shapes(obj['world']),
           len(obj['world']),
           "unclamped",
           "clamped")

    print(c_code)
    
if __name__ == '__main__':
    import sys
    x = yaml_file_to_world_objects(sys.argv[1])
    world_objects_to_c_file(x)


