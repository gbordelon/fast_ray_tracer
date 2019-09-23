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

        return """Material material_{0} = material_alloc();
    material_{0}->color[0] = {2:.10f};
    material_{0}->color[1] = {3:.10f};
    material_{0}->color[2] = {4:.10f};
    material_{0}->ambient = {5:.10f};
    material_{0}->diffuse = {6:.10f};
    material_{0}->specular = {7:.10f};
    material_{0}->shininess = {8:.10f};
    material_{0}->reflective = {9:.10f};
    material_{0}->transparency = {10:.10f};
    material_{0}->refractive_index = {11:.10f};
    material_{0}->casts_shadow = {12};
//    material_{0}->pattern = pattern_{0};
""".format(name,
           "",
           self.yaml_obj['color'][0], self.yaml_obj['color'][1], self.yaml_obj['color'][2],
           self.yaml_obj['ambient'],
           self.yaml_obj['diffuse'],
           self.yaml_obj['specular'],
           self.yaml_obj['shininess'],
           self.yaml_obj['reflective'],
           self.yaml_obj['transparency'],
           self.yaml_obj['refractive-index'],
           shadow_str)
