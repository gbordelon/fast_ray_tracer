from transform import Transform
from material import Material
from pattern import Pattern

from copy import deepcopy


class Shape(object):
    def __init__(self, yaml_obj, defines=None):
        self.yaml_obj = yaml_obj
        self.defines = defines
        self.name = ...

    @classmethod
    def from_yaml(cls, obj, defines) -> 'Shape':
        if 'transform' not in obj:
            obj['transform'] = []

        if obj["add"] == "sphere":
            return Sphere.from_yaml(obj)
        elif obj["add"] == "plane":
            return Plane.from_yaml(obj)
        elif obj["add"] == "cube":
            return Cube.from_yaml(obj)
        elif obj["add"] == "cone":
            return Cone.from_yaml(obj)
        elif obj["add"] == "cylinder":
            return Cylinder.from_yaml(obj)
        elif obj["add"] == "triangle":
            return Triangle.from_yaml(obj)
        elif obj["add"] == "smooth-triangle":
            return SmoothTriangle.from_yaml(obj)
        elif obj['add'] in ['toroid', 'torus']:
            return Toroid.from_yaml(obj);
        elif obj["add"] == "group":
            if "material" in obj:
                for child in obj["children"]:
                    if "material" not in child:
                        child["material"] = deepcopy(obj["material"])
            return Group.from_yaml(obj, defines)
        elif obj["add"] == "csg":
            if "material" in obj:
                if "material" not in obj["left"]:
                    obj["left"]["material"] = deepcopy(obj["material"])
                if "material" not in obj["right"]:
                    obj["right"]["material"] = deepcopy(obj["material"])
            return CSG.from_yaml(obj, defines)
        elif obj["add"] == "obj":
            from obj_parser import OBJParser
            return OBJParser.from_yaml(obj, defines)

        raise ValueError("unsupported shape.")

    def get_file_name(self):
        if 'file' in self.yaml_obj:
            return self.yaml_obj['file']
        return None

    def c_repr(self, name, parent_name, offset, resources):
        raise NotImplementedError

class CSG(Shape):
    def __init__(self, yaml_obj, defines):
        Shape.__init__(self, yaml_obj, defines)

    @classmethod
    def from_yaml(cls, tree, defines) -> 'CSG':
        return cls(tree, defines)

    def c_repr(self, name, parent_name, offset, resources):
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        left = Shape.from_yaml(self.yaml_obj['left'], self.defines)
        right = Shape.from_yaml(self.yaml_obj['right'], self.defines)

        if 'op' in self.yaml_obj:
            op = self.yaml_obj['op']
        else:
            op = self.yaml_obj['operation']

        if op == 'difference':
            op_str = 'CSG_DIFFERENCE'
        elif op == 'intersection':
            op_str = 'CSG_INTERSECT'
        elif op == 'union':
            op_str = 'CSG_UNION'
        else:
            raise ValueError("Unknown CSG operation: {}".format(op))

        buf = """
    /* children for {0} */
    Shape shape_{0}_children = array_of_shapes({1});
""".format(name, 2)

        buf += """    {0}""".format(left.c_repr('{}_left'.format(name), 'shape_{0}_children'.format(name), 0, resources))

        file_name = left.get_file_name()
        if file_name is not None:
            resources[file_name] = 'shape_{}_left'.format(name)

        buf += """    {0}""".format(right.c_repr('{}_right'.format(name), 'shape_{0}_children'.format(name), 1, resources))

        file_name = right.get_file_name()
        if file_name is not None:
            resources[file_name] = 'shape_{}_right'.format(name)


        buf += """
    /* end children for {0} */

    {2}
    Shape shape_{0} = {1} + {4};
    csg(shape_{0}, {3}, shape_{0}_left, shape_{0}_right);
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           parent_name,
           transform.c_repr(name),
           op_str,
           offset)

        return buf


class Group(Shape):
    def __init__(self, yaml_obj, defines):
        Shape.__init__(self, yaml_obj, defines)

    @classmethod
    def from_yaml(cls, tree, defines) -> 'Group':
        return cls(tree, defines)

    def c_repr(self, name, parent_name, offset, resources):
        if 'transform' not in self.yaml_obj:
            yaml_obj['transform'] = []

        transform = Transform.from_yaml(self.yaml_obj['transform'])
        children = []
        for child in self.yaml_obj['children']:
            children.append(Shape.from_yaml(child, self.defines))
            
        buf = """
    /* children for {0} */
    Shape shape_{0}_children = array_of_shapes({1});
""".format(name, len(children))

        for i, child in enumerate(children):
            buf += """
    {0}
""".format(child.c_repr('{0}_child_{1}'.format(name, i), 'shape_{0}_children'.format(name), i, resources))

            file_name = child.get_file_name()
            if file_name is not None:
                resources[file_name] = 'shape_{0}_child_{1}'.format(name, i)
        buf += """    /* end children for {0} */

    {2}
    Shape shape_{0} = {1} + {4};
    group(shape_{0}, shape_{0}_children, {3});
    //shape_free(shape_{0}_children);
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
          parent_name,
          transform.c_repr(name),
          len(children),
          offset)

        return buf


class Sphere(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'Sphere':
        return cls(obj)

    def c_repr(self, name, parent_name, offset, resources):
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {2}
    {3}
    Shape shape_{0} = {1} + {4};
    sphere(shape_{0});
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           parent_name,
           material.c_repr(name),
           transform.c_repr(name),
           offset)

        return buf

class Toroid(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'Toroid':
        if 'r1' not in obj:
            obj['r1'] = 0.75
        if 'r2' not in obj:
            obj['r2'] = 0.25
        return cls(obj)

    def c_repr(self, name, parent_name, offset, resources):
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {2}
    {3}
    Shape shape_{0} = {1} + {4};
    toroid(shape_{0});
    shape_{0}->fields.toroid.r1 = {5};
    shape_{0}->fields.toroid.r2 = {6};
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           parent_name,
           material.c_repr(name),
           transform.c_repr(name),
           offset,
           self.yaml_obj['r1'],
           self.yaml_obj['r2'])

        return buf

class Plane(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'Plane':
        return cls(yaml_obj=obj)

    def c_repr(self, name, parent_name, offset, resources):
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {2}
    {3}
    Shape shape_{0} = {1} + {4};
    plane(shape_{0});
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           parent_name,
           material.c_repr(name),
           transform.c_repr(name),
           offset)

        return buf


class Cube(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'Cube':
        return cls(yaml_obj=obj)

    def c_repr(self, name, parent_name, offset, resources):
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {2}
    {3}
    Shape shape_{0} = {1} + {4};
    cube(shape_{0});
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           parent_name,
           material.c_repr(name),
           transform.c_repr(name),
           offset)

        return buf


class Cone(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'Cone':
        return cls(yaml_obj=obj)

    def c_repr(self, name, parent_name, offset, resources):
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {2}
    {3}
    Shape shape_{0} = {1} + {4};
    cone(shape_{0});
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           parent_name,
           material.c_repr(name),
           transform.c_repr(name),
           offset)

        if 'min' in self.yaml_obj:
            buf += "    shape_{0}->fields.cone.minimum = {1:.10f};\n".format(name, self.yaml_obj['min'])
        if 'max' in self.yaml_obj:
            buf += "    shape_{0}->fields.cone.maximum = {1:.10f};\n".format(name, self.yaml_obj['max'])
        if 'closed' in self.yaml_obj:
            bool_str = 'false'
            if self.yaml_obj['closed']:
                bool_str = 'true'
            buf += "    shape_{0}->fields.cone.closed = {1};\n".format(name, bool_str)

        return buf


class Cylinder(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'Cylinder':
        b = True
        if 'closed' in obj:
            b = obj['closed']
        return cls(yaml_obj=obj)

    def c_repr(self, name, parent_name, offset, resources):
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {2}
    {3}
    Shape shape_{0} = {1} + {4};
    cylinder(shape_{0});
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});
""".format(name,
           parent_name,
           material.c_repr(name),
           transform.c_repr(name),
           offset)

        if 'min' in self.yaml_obj:
            buf += "    shape_{0}->fields.cylinder.minimum = {1:.10f};\n".format(name, self.yaml_obj['min'])
        if 'max' in self.yaml_obj:
            buf += "    shape_{0}->fields.cylinder.maximum = {1:.10f};\n".format(name, self.yaml_obj['max'])
        if 'closed' in self.yaml_obj:
            bool_str = 'false'
            if self.yaml_obj['closed']:
                bool_str = 'true'
            buf += "    shape_{0}->fields.cylinder.closed = {1};\n".format(name, bool_str)

        return buf


class Triangle(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'Triangle':
        return cls(yaml_obj=obj)

    def c_repr(self, name, parent_name, offset, resources):
        if 'material' not in self.yaml_obj:
            self.yaml_obj['material'] = {}
        if 'transform' not in self.yaml_obj:
            self.yaml_obj['transform'] = []
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {1}
    {2}
    Point shape_{0}_p1 = point({3:.10f}, {4:.10f}, {5:.10f});
    Point shape_{0}_p2 = point({6:.10f}, {7:.10f}, {8:.10f});
    Point shape_{0}_p3 = point({9:.10f}, {10:.10f}, {11:.10f});

    Shape shape_{0} = {12} + {13};
    triangle(shape_{0}, shape_{0}_p1->arr, shape_{0}_p2->arr, shape_{0}_p3->arr);
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});

    point_free(shape_{0}_p3);
    point_free(shape_{0}_p2);
    point_free(shape_{0}_p1);
""".format(name,
          material.c_repr(name),
          transform.c_repr(name),
          self.yaml_obj['p1'][0], self.yaml_obj['p1'][1], self.yaml_obj['p1'][2],
          self.yaml_obj['p2'][0], self.yaml_obj['p2'][1], self.yaml_obj['p2'][2],
          self.yaml_obj['p3'][0], self.yaml_obj['p3'][1], self.yaml_obj['p3'][2],
          parent_name, offset)

        return buf

class SmoothTriangle(Shape):
    def __init__(self, yaml_obj):
        Shape.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'SmoothTriangle':
        return cls(yaml_obj=obj)

    def c_repr(self, name, parent_name, offset, resources):
        if 'material' not in self.yaml_obj:
            self.yaml_obj['material'] = {}
        if 'transform' not in self.yaml_obj:
            self.yaml_obj['transform'] = []
        material = Material.from_yaml(self.yaml_obj['material'])
        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = """
    {1}
    {2}
    Point shape_{0}_p1 = point({3:.10f}, {4:.10f}, {5:.10f});
    Point shape_{0}_p2 = point({6:.10f}, {7:.10f}, {8:.10f});
    Point shape_{0}_p3 = point({9:.10f}, {10:.10f}, {11:.10f});
    Vector shape_{0}_n1 = vector({12:.10f}, {13:.10f}, {14:.10f});
    Vector shape_{0}_n2 = vector({15:.10f}, {16:.10f}, {17:.10f});
    Vector shape_{0}_n3 = vector({18:.10f}, {19:.10f}, {20:.10f});

    Shape shape_{0} = {21} + {22};
    smooth_triangle(shape_{0}, shape_{0}_p1->arr, shape_{0}_p2->arr, shape_{0}_p3->arr, shape_{0}_n1->arr, shape_{0}_n2->arr, shape_{0}_n3->arr);
    shape_set_material(shape_{0}, material_{0});
    shape_set_transform(shape_{0}, transform_{0});

    vector_free(shape_{0}_n3);
    vector_free(shape_{0}_n2);
    vector_free(shape_{0}_n1);
    point_free(shape_{0}_p3);
    point_free(shape_{0}_p2);
    point_free(shape_{0}_p1);
    
""".format(name,
          material.c_repr(name),
          transform.c_repr(name),
          self.yaml_obj['p1'][0], self.yaml_obj['p1'][1], self.yaml_obj['p1'][2],
          self.yaml_obj['p2'][0], self.yaml_obj['p2'][1], self.yaml_obj['p2'][2],
          self.yaml_obj['p3'][0], self.yaml_obj['p3'][1], self.yaml_obj['p3'][2],
          self.yaml_obj['n1'][0], self.yaml_obj['n1'][1], self.yaml_obj['n1'][2],
          self.yaml_obj['n2'][0], self.yaml_obj['n2'][1], self.yaml_obj['n2'][2],
          self.yaml_obj['n3'][0], self.yaml_obj['n3'][1], self.yaml_obj['n3'][2],
          parent_name, offset)

        return buf


def allocate_shapes(list_of_shapes):
    num = len(list_of_shapes)
    resources = {}
    buf = """    /* shapes */
    Shape all_shapes = array_of_shapes({0});
""".format(num)
    for i, shape in enumerate(list_of_shapes):
        buf += """
    /* shape {0} */
{1}
    /* end shape {0} */""".format(i, shape.c_repr(i, "all_shapes", i, resources))
    file_name = shape.get_file_name()
    if file_name is not None:
        resources[file_name] = 'shape_{0}'.format(i)

    buf += """
    /* end shapes */
"""

    return buf
