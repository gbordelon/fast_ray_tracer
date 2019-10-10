from pattern import Pattern

class Material(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Material':
        c = [1,1,1]
        ambient = 0.1
        diffuse = 0.9
        specular = 0.9
        shininess = 200.0
        reflective = 0.0
        transparency = 0.0
        refractive_index = 1.0
        casts_shadow = True
        pat = None

        if obj is None:
            obj = {}

        if 'color' not in obj:
            obj['color'] = c

        if 'diffuse' not in obj:
            obj['diffuse'] = diffuse

        if 'ambient' not in obj:
            obj['ambient'] = ambient

        if 'specular' not in obj:
            obj['specular'] = specular

        if 'shininess' not in obj:
            obj['shininess'] = shininess

        if 'reflective' not in obj:
            obj['reflective'] = reflective

        if 'transparency' not in obj:
            obj['transparency'] = transparency

        if 'refractive-index' not in obj:
            obj['refractive-index'] = refractive_index

        if 'shadow' not in obj:
            obj['shadow'] = casts_shadow


        return Material(yaml_obj=obj)

    def c_repr(self, name):
        shadow_str = 'false'
        if self.yaml_obj['shadow']:
            shadow_str = 'true'

        if 'pattern' not in self.yaml_obj:
            pat = Pattern.from_yaml({})
        else:
            pat = Pattern.from_yaml(self.yaml_obj['pattern'])
        return """{1}
    Color material_{0}_color_raw = color({2:.10f}, {3:.10f}, {4:.10f});

    Material material_{0} = material_alloc();
    color_space_fn(material_{0}_color_raw, material_{0}->color);

    material_{0}->ambient = {5:.10f};
    material_{0}->diffuse = {6:.10f};
    material_{0}->specular = {7:.10f};
    material_{0}->shininess = {8:.10f};
    material_{0}->reflective = {9:.10f};
    material_{0}->transparency = {10:.10f};
    material_{0}->refractive_index = {11:.10f};
    material_{0}->casts_shadow = {12};
    material_{0}->pattern = pattern_{0};
""".format(name,
           pat.c_repr(name),
           self.yaml_obj['color'][0], self.yaml_obj['color'][1], self.yaml_obj['color'][2],
           self.yaml_obj['ambient'],
           self.yaml_obj['diffuse'],
           self.yaml_obj['specular'],
           self.yaml_obj['shininess'],
           self.yaml_obj['reflective'],
           self.yaml_obj['transparency'],
           self.yaml_obj['refractive-index'],
           shadow_str)
