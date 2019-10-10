class GlobalConfig(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        if 'illumination' not in obj:
            obj['illumination'] = {}
        if 'threading' not in obj:
            obj['threading'] = {}
        if 'scene' not in obj:
            obj['scene'] = {}
        return cls(yaml_obj=obj)

    def c_repr(self):
        illum = IlluminationConfig.from_yaml(self.yaml_obj['illumination'])
        threading = ThreadingConfig.from_yaml(self.yaml_obj['threading'])
        scene = SceneConfig.from_yaml(self.yaml_obj['scene'])
        return """    /* config */
    struct global_config global_config;
{0}
{1}
{2}
    /* end config */
""".format(illum.c_repr(),
           threading.c_repr(),
           scene.c_repr())

class SceneConfig(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        if 'divide-threshold' not in obj:
            obj['divide-threshold'] = 1
        return cls(yaml_obj=obj)

    def c_repr(self):
        return """    global_config.scene.divide_threshold = {0};""".format(self.yaml_obj['divide-threshold'])

class ThreadingConfig(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        if 'thread-count' not in obj:
            obj['thread-count'] = 4
        return cls(yaml_obj=obj)

    def c_repr(self):
        return """    global_config.threading.num_threads = {0};""".format(self.yaml_obj['thread-count'])

class DirectIlluminationConfig(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        if 'include-ambient' not in obj:
            obj['include-ambient'] = True
        if 'include-diffuse' not in obj:
            obj['include-diffuse'] = True
        if 'include-specular-highlight' not in obj:
            obj['include-specular-highlight'] = True
        if 'include-specular' not in obj:
            obj['include-specular'] = True

        return cls(yaml_obj=obj)

    def c_repr(self):
        return """    global_config.illumination.di.include_ambient = {0};
    global_config.illumination.di.include_diffuse = {1};
    global_config.illumination.di.include_specular_highlight = {2};
    global_config.illumination.di.include_specular = {3};
""".format(self.yaml_obj['include-ambient'],
           self.yaml_obj['include-diffuse'],
           self.yaml_obj['include-specular-highlight'],
           self.yaml_obj['include-specular'])

class GlobalIlluminationConfig(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        if 'include-caustics' not in obj:
            obj['include-caustics'] = False
        if 'include-final-gather' not in obj:
            obj['include-final-gather'] = False
        if 'usteps' not in obj:
            obj['usteps'] = 1
        if 'vsteps' not in obj:
            obj['vsteps'] = 1
        if 'irradiance-estimate-num' not in obj:
            obj['irradiance-estimate-num'] = 200
        if 'irradiance-estimate-radius' not in obj:
            obj['irradiance-estimate-radius'] = 0.1
        if 'photon-count' not in obj:
            obj['photon-count'] = 0

        return cls(yaml_obj=obj)

    def c_repr(self):
        return """    global_config.illumination.gi.include_caustics = {0};
    global_config.illumination.gi.include_final_gather = {1};
    global_config.illumination.gi.usteps = {2};
    global_config.illumination.gi.vsteps = {3};
    global_config.illumination.gi.irradiance_estimate_num = {4};
    global_config.illumination.gi.irradiance_estimate_radius = {5:.10f};
    global_config.illumination.gi.photon_count = {6};
""".format(self.yaml_obj['include-caustics'],
           self.yaml_obj['include-final-gather'],
           self.yaml_obj['usteps'],
           self.yaml_obj['vsteps'],
           self.yaml_obj['irradiance-estimate-num'],
           self.yaml_obj['irradiance-estimate-radius'],
           self.yaml_obj['photon-count'])

class IlluminationConfig(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        if 'include-direct' not in obj:
            obj['include-direct'] = True
        if 'include-global' not in obj:
            obj['include-global'] = False
        if 'visualize-photon-map' not in obj:
            obj['visualize-photon-map'] = False
        if 'visualize-soft-indirect' not in obj:
            obj['visualize-soft-indirect'] = False
        if 'global-illumination' not in obj:
            obj['global-illumination'] = {}
        if 'direct-illumination' not in obj:
            obj['direct-illumination'] = {}

        return cls(yaml_obj=obj)

    def c_repr(self):
        gi = GlobalIlluminationConfig.from_yaml(self.yaml_obj['global-illumination'])
        di = DirectIlluminationConfig.from_yaml(self.yaml_obj['direct-illumination'])
        return """    global_config.illumination.include_direct = {0};
    global_config.illumination.include_global = {1};
    global_config.illumination.debug_visualize_photon_map = {2};
    global_config.illumination.debug_visualize_soft_indirect = {3};
{4}
{5}""".format(self.yaml_obj['include-direct'],
           self.yaml_obj['include-global'],
           self.yaml_obj['visualize-photon-map'],
           self.yaml_obj['visualize-soft-indirect'],
           di.c_repr(),
           gi.c_repr())

