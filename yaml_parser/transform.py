class Transform(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        return cls(yaml_obj=obj)

    @classmethod
    def single_transform(cls, item, name):
        if item[0] == 'translate':
            return "matrix_translate({:.10f}, {:.10f}, {:.10f}, {})".format(item[1], item[2], item[3], name)
        elif item[0] == 'scale':
            return "matrix_scale({:.10f}, {:.10f}, {:.10f}, {})".format(item[1], item[2], item[3], name)
        elif item[0] == 'rotate-x':
            return "matrix_rotate_x({:.10f}, {})".format(item[1], name)
        elif item[0] == 'rotate-y':
            return "matrix_rotate_y({:.10f}, {})".format(item[1], name)
        elif item[0] == 'rotate-z':
            return "matrix_rotate_z({:.10f}, {})".format(item[1], name)
        elif item[0] == 'shear':
            return "matrix_shear({:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f}, {})".format(item[1], item[2], item[3], item[4], item[5], item[6], name)
        else:
            raise ValueError('Unknown translation type: {}'.format(item[0]))

    def c_repr(self, name):
        buf = "Matrix transform_{0}".format(name)
        if len(self.yaml_obj) == 0:
            buf += ";\n    matrix_identity(transform_{0});".format(name)
        elif len(self.yaml_obj) > 1:
            buf += """, transform_{0}_tmp;
    matrix_identity(transform_{0});
""".format(name)
            for m in self.yaml_obj:
                buf += """    {1};
    transform_chain(transform_{0}_tmp, transform_{0});
""".format(name, Transform.single_transform(m, "transform_{0}_tmp".format(name)))
        else:
            buf += ";\n    {};".format(Transform.single_transform(self.yaml_obj[-1], "transform_{0}".format(name)))
        return buf
