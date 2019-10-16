from shapes import Material, Shape
from transform import Transform

class OBJParser(Shape):
    def __init__(self, yaml_obj, defines=None):
        Shape.__init__(self, yaml_obj, defines)
        self.yaml_obj = yaml_obj
        self.defines = defines

    @classmethod
    def from_yaml(cls, obj, defines) -> 'OBJParser':
        return cls(obj, defines)

    def c_repr(self, name, parent_name, offset, resources):
        material = None
        if 'transform' not in self.yaml_obj:
            self.yaml_obj['transform'] = []
        if 'material' in self.yaml_obj:
            material = Material.from_yaml(self.yaml_obj['material'])

        transform = Transform.from_yaml(self.yaml_obj['transform'])
        file_path = self.yaml_obj['file']

        buf = ""
        if material is not None:
            buf += "{0}".format( material.c_repr(name, resources))
        buf += """
    {3}
    Shape shape_{0} = {1} + {2};
""".format(name, parent_name, offset, transform.c_repr(name))
        if file_path in resources:
            buf += """    shape_copy({2}, NULL, shape_{0});""".format(name, file_path, resources[file_path])
        else:
            buf += """
    if (access("{1}", F_OK ) == -1 ) {{
        printf("file '{1}' does not exist.");
        return 1;
    }}
    printf("Loading resource '{1}'... ");
    fflush(stdout);
    construct_group_from_obj_file("{1}", color_space_fn, shape_{0});
    printf("Done!\\n");
    fflush(stdout);
""".format(name, file_path, file_path[:-3] + 'mtl')

        if material is not None:
            buf += """
    shape_set_material_recursive(shape_{0}, material_{0});""".format(name)
        buf += """
    shape_set_transform(shape_{0}, transform_{0});
""".format(name)

        return buf

def obj_to_group(parser):
    return obj_group

