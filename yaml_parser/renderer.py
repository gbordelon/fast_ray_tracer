class Aperture(object):
    def __init__(self, yaml_obj=None):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Aperture':
        if 'size' not in obj:
            obj['size'] = 0.0
        if 'type' not in obj:
            obj['type'] = [ 'POINT_APERTURE' ]
        if 'jitter' not in obj:
            obj['jitter'] = False
        return cls(yaml_obj=obj)

    def c_repr(self):
        bool_str = "false"
        if self.yaml_obj['jitter']:
            bool_str = "true"

        buf = """    struct aperture aperture;
    aperture.type = {0};
    aperture.size = {1};
    aperture.jitter = {2};
""".format(self.yaml_obj['type'][0],
           self.yaml_obj['size'],
           bool_str)

        if self.yaml_obj['type'][0] == 'CIRCULAR_APERTURE':
            buf += """    aperture.u.circle.r1 = {0:.10f};
""".format(*self.yaml_obj['type'][1:])
        elif self.yaml_obj['type'][0] == 'CROSS_APERTURE':
            buf += """    aperture.u.cross.x1 = {0:.10f};
    aperture.u.cross.x2 = {1:.10f};
    aperture.u.cross.y1 = {2:.10f};
    aperture.u.cross.y2 = {3:.10f};
""".format(*self.yaml_obj['type'][1:])
        elif self.yaml_obj['type'][0] == 'DIAMOND_APERTURE':
            buf += """    aperture.u.diamond.b1 = {0:.10f};
    aperture.u.diamond.b2 = {1:.10f};
    aperture.u.diamond.b3 = {2:.10f};
    aperture.u.diamond.b4 = {3:.10f};
""".format(*self.yaml_obj['type'][1:])
        elif self.yaml_obj['type'][0] == 'DOUBLE_CIRCLE_APERTURE':
            buf += """    aperture.u.double_circle.r1 = {0:.10f};
    aperture.u.double_circle.r2 = {1:.10f};
""".format(*self.yaml_obj['type'][1:])
        else:
            pass

        return buf

class Camera(object):
    def __init__(self, yaml_obj=None):
       self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Camera':
        if 'focal-length' not in obj:
            obj['focal-length'] = 1.0
        if 'samples-per-pixel' not in obj:
            obj['samples-per-pixel'] = 1
        if 'aperture' not in obj:
            obj['aperture'] = {}
        return cls(yaml_obj=obj)

    def c_repr(self):
        aperture = Aperture.from_yaml(self.yaml_obj['aperture'])
        return """    /* camera */
{0}
    Point from = point({1:.10f}, {2:.10f}, {3:.10f});
    Point to = point({4:.10f}, {5:.10f}, {6:.10f});
    Vector up = vector({7:.10f}, {8:.10f}, {9:.10f});

    Camera cam = camera({10}, {11}, {12:.10f}/*field_of_view*/, {13:.10f}/*distance*/, aperture, {14}/*sample_num*/, view_transform(from, to, up));

    vector_free(up);
    point_free(to);
    point_free(from);
    /* end camera */
""".format(aperture.c_repr(),
           self.yaml_obj['from'][0], self.yaml_obj['from'][1], self.yaml_obj['from'][2],
           self.yaml_obj['to'][0], self.yaml_obj['to'][1], self.yaml_obj['to'][2],
           self.yaml_obj['up'][0], self.yaml_obj['up'][1], self.yaml_obj['up'][2],
           self.yaml_obj['width'], self.yaml_obj['height'], self.yaml_obj['field-of-view'],
           self.yaml_obj['focal-length'], self.yaml_obj['samples-per-pixel'])

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

    def c_repr(self, name):
        return """    /* point light {0} */
    Light point_light_{0} = all_lights + {0};
    Point point_light_{0}_point = point({1:.10f}, {2:.10f}, {3:.10f});
    Color point_light_{0}_intensity = color({4:.10f}, {5:.10f}, {6:.10f});
    point_light(point_light_{0}_point, point_light_{0}_intensity, point_light_{0});

    color_free(point_light_{0}_intensity);
    point_free(point_light_{0}_point);
    /* end point light {0} */
""".format(name,
           self.yaml_obj['at'][0], self.yaml_obj['at'][1], self.yaml_obj['at'][2],
           self.yaml_obj['intensity'][0], self.yaml_obj['intensity'][1], self.yaml_obj['intensity'][2])

class AreaLight(Light):
    def __init__(self, yaml_obj):
        Light.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'AreaLight':
        return cls(obj)

    def c_repr(self, name):
        bool_str = "false"
        if self.yaml_obj['jitter']:
            bool_str = "true"
        return """
    /* area light {0} */
    Light area_light_{0} = all_lights + {0};
    Point area_light_{0}_corner = point({1:.10f}, {2:.10f}, {3:.10f});
    Color area_light_{0}_intensity = color({4:.10f}, {5:.10f}, {6:.10f});
    Vector area_light_{0}_uvec = vector({7:.10f}, {8:.10f}, {9:.10f});
    Vector area_light_{0}_vvec = vector({10:.10f}, {11:.10f}, {12:.10f});
    area_light(area_light_{0}_corner->arr, area_light_{0}_uvec->arr, {13}/*usteps*/, area_light_{0}_vvec->arr, {14}/*vsteps*/, {15}/*jitter*/, area_light_{0}_intensity->arr, area_light_{0});

    vector_free(area_light_{0}_vvec);
    vector_free(area_light_{0}_uvec);
    color_free(area_light_{0}_intensity);
    point_free(area_light_{0}_corner);
    /* end area light 0 */
""".format(name,
           self.yaml_obj['corner'][0], self.yaml_obj['corner'][1], self.yaml_obj['corner'][2],
           self.yaml_obj['intensity'][0], self.yaml_obj['intensity'][1], self.yaml_obj['intensity'][2],
           self.yaml_obj['uvec'][0], self.yaml_obj['uvec'][1], self.yaml_obj['uvec'][2],
           self.yaml_obj['vvec'][0], self.yaml_obj['vvec'][1], self.yaml_obj['vvec'][2],
           self.yaml_obj['usteps'], self.yaml_obj['vsteps'], bool_str)



def allocate_lights(list_of_lights):
    num = len(list_of_lights)
    buf = """    /* lights */
    Light all_lights = array_of_lights({0});

""".format(num);
    for i, light in enumerate(list_of_lights):
        buf += light.c_repr(i)

    buf += """
    /* end lights */
"""
    return buf;
