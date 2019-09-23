class Transform(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        return cls(yaml_obj=obj)

    @classmethod
    def single_transform(cls, item):
        if item[0] == 'translate':
            return "matrix_translate_alloc({:.10f}, {:.10f}, {:.10f})".format(item[1], item[2], item[3])
        elif item[0] == 'scale':
            return "matrix_scale_alloc({:.10f}, {:.10f}, {:.10f})".format(item[1], item[2], item[3])
        elif item[0] == 'rotate-x':
            return "matrix_rotate_x_alloc({:.10f})".format(item[1])
        elif item[0] == 'rotate-y':
            return "matrix_rotate_y_alloc({:.10f})".format(item[1])
        elif item[0] == 'rotate-z':
            return "matrix_rotate_z_alloc({:.10f})".format(item[1])
        elif item[0] == 'shear':
            return "matrix_shear_alloc({:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f}, {:.10f})".format(item[1], item[2], item[3], item[4], item[5], item[6])
        else:
            raise ValueError('Unknown translation type: {}'.format(item[0]))



##
## should be
##    transform_chain(transform_chain(matrix_rotate_y_alloc(), matrix_rotate_x_alloc()),
##                    matrix_scale_alloc());
##
## was
##    ransform_chain(transform_chain(matrix_rotate_y_alloc(), matrix_scale_alloc()),
##                   matrix_rotate_x_alloc());
##  
    def c_repr(self, name):
        buf = "Matrix transform_{0} = ".format(name)
        if len(self.yaml_obj) == 0:
            buf += "matrix_identity_alloc();"
        elif len(self.yaml_obj) > 1:
            for m in self.yaml_obj[:-1]:
                buf += "transform_chain("
            buf += Transform.single_transform(self.yaml_obj[-1])
            for m in reversed(self.yaml_obj[:-1]):
                buf += ", {})".format(Transform.single_transform(m))
            buf += ";"
        else:
            buf += "{};".format(Transform.single_transform(self.yaml_obj[-1]))
        return buf
