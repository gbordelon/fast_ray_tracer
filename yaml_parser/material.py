from pattern import Pattern

map_types = ['Ka', 'Kd', 'Ks', 'Ns', 'bump', 'disp']

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

    def c_repr(self, name, resources):
        shadow_str = 'false'
        if self.yaml_obj['shadow']:
            shadow_str = 'true'

        if 'pattern' not in self.yaml_obj and 'patterns' not in self.yaml_obj:
            pats = {x : Pattern.from_yaml({}) for x in map_types}
        elif 'pattern' in self.yaml_obj:
            pats = {}
            pats = {x : Pattern.from_yaml({}) for x in map_types}
            pats['Ka'] = Pattern.from_yaml(self.yaml_obj['pattern'])
            pats['Kd'] = pats['Ka']
        elif 'patterns' in self.yaml_obj:
            pats = {}
            pats = {x : Pattern.from_yaml({}) for x in map_types}
            for k in map_types:
                if k in self.yaml_obj['patterns']:
                    pats[k] = Pattern.from_yaml(self.yaml_obj['patterns'][k])
                
        buf = ""
        for pat in pats:
            buf += """{}""".format(pats[pat].c_repr("{0}_{1}".format(name, pat), resources))
        return buf + """    Color material_{0}_color_raw = color({2:.10f}, {3:.10f}, {4:.10f});

    Material material_{0} = material_alloc();
    color_space_fn(material_{0}_color_raw, material_{0}->Ka);
    color_space_fn(material_{0}_color_raw, material_{0}->Kd);
    color_space_fn(material_{0}_color_raw, material_{0}->Ks);
    color_scale(material_{0}->Ka, {5:.10f});
    color_scale(material_{0}->Kd, {6:.10f});
    color_scale(material_{0}->Ks, {7:.10f});
    material_{0}->reflective = {9:.10f};

    material_{0}->Ns = {8:.10f};
    material_{0}->Tr = {10:.10f};
    material_{0}->Tf[0] = {10:.10f};
    material_{0}->Tf[1] = {10:.10f};
    material_{0}->Tf[2] = {10:.10f};
    material_{0}->Ni = {11:.10f};
    material_{0}->casts_shadow = {12};
    material_set_pattern(material_{0}, map_Ka, pattern_{0}_Ka);
    material_set_pattern(material_{0}, map_Kd, pattern_{0}_Kd);
    material_set_pattern(material_{0}, map_Ks, pattern_{0}_Ks);
    material_set_pattern(material_{0}, map_Ns, pattern_{0}_Ns);
    material_set_pattern(material_{0}, map_bump, pattern_{0}_bump);
    material_set_pattern(material_{0}, map_disp, pattern_{0}_disp);
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
