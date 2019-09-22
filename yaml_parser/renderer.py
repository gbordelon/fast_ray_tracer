class Camera(object):
    def __init__(self, yaml_obj=None):
       self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Camera':
        if 'aperture-size' not in obj:
            obj['aperture-size'] = 0.0
        if 'focal-length' not in obj:
            obj['focal-length'] = 1.0
        if 'aperture-shape' not in obj:
            obj['aperture-shape'] = 'POINT_APERTURE'
        if 'samples-per-pixel' not in obj:
            obj['samples-per-pixel'] = 1
        return cls(yaml_obj=obj)

    def __repr__(self):
        return """
    /* camera */
    Point from = point({:.10f}, {:.10f}, {:.10f});
    Point to = point({:f}, {:f}, {:f});
    Vector up = vector({:f}, {:f}, {:f});
    Camera cam = camera({}, {}, {:f}/*field_of_view*/, {:f}/*aperture*/, {:f}/*distance*/, {}/*aperture shape*/, {}/*sample_num*/, view_transform(from, to, up));
    /* end camera */
""".format(self.yaml_obj['from'][0], self.yaml_obj['from'][1], self.yaml_obj['from'][2],
           self.yaml_obj['to'][0], self.yaml_obj['to'][1], self.yaml_obj['to'][2],
           self.yaml_obj['up'][0], self.yaml_obj['up'][1], self.yaml_obj['up'][2],
           self.yaml_obj['width'], self.yaml_obj['height'], self.yaml_obj['field-of-view'], self.yaml_obj['aperture-size'], self.yaml_obj['focal-length'], self.yaml_obj['aperture-shape'], self.yaml_obj['samples-per-pixel'])

class Light(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Light':
        if 'at' in obj:
            return PointLight.from_yaml(obj)
        elif 'corner' in obj:
            return AreaLight.from_yaml(obj)

class PointLight(Light):
    def __init__(self, yaml_obj):
        Light.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'PointLight':
        return cls(obj)

    def __repr__(self):
        return """
    /* point light 0 */
    Point point_light_0_point = point({:.10f}, {:.10f}, {:.10f});
    Color point_light_0_intensity = color({:.10f}, {:.10f}, {:.10f});
    Light point_light_0 = point_light_alloc(point_light_0, point_light_0_intensity);
    /* end point light 0 */
""".format(self.yaml_obj['at'][0], self.yaml_obj['at'][1], self.yaml_obj['at'][2],
           self.yaml_obj['intensity'][0], self.yaml_obj['intensity'][1], self.yaml_obj['intensity'][2])

class AreaLight(Light):
    def __init__(self, yaml_obj):
        Light.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'AreaLight':
        return cls(obj)

    def __repr__(self):
        return """
    /* area light 0 */
    Point area_light_0_corner = point({:.10f}, {:.10f}, {:.10f});
    Color area_light_0_intensity = color({:.10f}, {:.10f}, {:.10f});
    Vector area_light_0_uvec = vector({:.10f}, {:.10f}, {:.10f});
    Vector area_light_0_vvec = vector({:.10f}, {:.10f}, {:.10f});
    Light area_light_0 = area_light_alloc(area_light_0_corner->arr, area_light_0_uvec->arr, {}/*usteps*/, area_light_0_vvec->arr, {}/*vsteps*/, {}/*jitter*/, area_light_0_intensity->arr);
    /* end area light 0 */
""".format(self.yaml_obj['corner'][0], self.yaml_obj['corner'][1], self.yaml_obj['corner'][2],
           self.yaml_obj['intensity'][0], self.yaml_obj['intensity'][1], self.yaml_obj['intensity'][2],
           self.yaml_obj['uvec'][0], self.yaml_obj['uvec'][1], self.yaml_obj['uvec'][2],
           self.yaml_obj['vvec'][0], self.yaml_obj['vvec'][1], self.yaml_obj['vvec'][2],
           self.yaml_obj['usteps'], self.yaml_obj['vsteps'], self.yaml_obj['jitter'])
