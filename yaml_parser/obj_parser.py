from shapes import Material
from transform import Transform

class OBJParser(object):
    def __init__(self, yaml_obj, defines=None):
        self.yaml_obj = yaml_obj
        self.defines = defines

    @classmethod
    def from_yaml(cls, obj, defines) -> 'OBJParser':
        return cls(obj, defines)

    def c_repr(self, name, parent_name, offset):
        if 'material' not in self.yaml_obj:
            self.yaml_obj['material'] = {}
        if 'transform' not in self.yaml_obj:
            self.yaml_obj['transform'] = []
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        file_path = self.yaml_obj['file']

        buf = """{6}
    {7}
    Shape shape_{0} = {4} + {5};
    
    if (access("{1}", F_OK ) == -1 ) {{
        printf("file '{1}' does not exist.");
        return 1;
    }}
    printf("Loading resource '{1}'... ");
    fflush(stdout);
    construct_group_from_obj_file("{1}", shape_{0});
    printf("Done!\\n");
    fflush(stdout);
    shape_set_material_recursive(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           file_path,
           material.c_repr(name),
           transform.c_repr(name),
           parent_name,
           offset,
           material.c_repr(name),
           transform.c_repr(name))


        return buf


def obj_to_group(parser):
    return obj_group

