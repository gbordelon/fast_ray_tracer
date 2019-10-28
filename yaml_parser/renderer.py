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

        buf = """    struct aperture ap;
    aperture({0}, {1}, {2}, {3}, {4}, &ap);
""".format(self.yaml_obj['type'][0],
           self.yaml_obj['size'],
           self.yaml_obj['usteps'],
           self.yaml_obj['vsteps'],
           bool_str)

        if self.yaml_obj['type'][0] == 'CIRCULAR_APERTURE':
            buf += """    ap.u.circle.r1 = {0:.10f};
""".format(*self.yaml_obj['type'][1:])
        elif self.yaml_obj['type'][0] == 'CROSS_APERTURE':
            buf += """    ap.u.cross.x1 = {0:.10f};
    ap.u.cross.x2 = {1:.10f};
    ap.u.cross.y1 = {2:.10f};
    ap.u.cross.y2 = {3:.10f};
""".format(*self.yaml_obj['type'][1:])
        elif self.yaml_obj['type'][0] == 'DIAMOND_APERTURE':
            buf += """    ap.u.diamond.b1 = {0:.10f};
    ap.u.diamond.b2 = {1:.10f};
    ap.u.diamond.b3 = {2:.10f};
    ap.u.diamond.b4 = {3:.10f};
""".format(*self.yaml_obj['type'][1:])
        elif self.yaml_obj['type'][0] == 'DOUGHNUT_APERTURE':
            buf += """    ap.u.doughnut.r1 = {0:.10f};
    ap.u.doughnut.r2 = {1:.10f};
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
        if 'aperture' not in obj:
            obj['aperture'] = {}
        if 'usteps' not in obj:
            obj['usteps'] = 1
        if 'vsteps' not in obj:
            obj['vsteps'] = 1
        return cls(yaml_obj=obj)

    def c_repr(self):
        if 'usteps' not in self.yaml_obj['aperture']:
            self.yaml_obj['aperture']['usteps'] = self.yaml_obj['usteps']
        if 'vsteps' not in self.yaml_obj['aperture']:
            self.yaml_obj['aperture']['vsteps'] = self.yaml_obj['vsteps']
        aperture = Aperture.from_yaml(self.yaml_obj['aperture'])
        return """    /* camera */
{0}
    Point from = {{ {1:.10f}, {2:.10f}, {3:.10f}, 1.0 }};
    Point to = {{ {4:.10f}, {5:.10f}, {6:.10f}, 1.0 }};
    Vector up = {{ {7:.10f}, {8:.10f}, {9:.10f}, 0.0 }};
    Matrix camera_xform;
    view_transform(from, to, up, camera_xform);

    Camera cam = camera({10}, {11}, {12:.10f}/*field_of_view*/, {13:.10f}/*distance*/, {14}/*usteps*/, {15}/*vsteps*/, &ap, camera_xform);

    /* end camera */
""".format(aperture.c_repr(),
           self.yaml_obj['from'][0], self.yaml_obj['from'][1], self.yaml_obj['from'][2],
           self.yaml_obj['to'][0], self.yaml_obj['to'][1], self.yaml_obj['to'][2],
           self.yaml_obj['up'][0], self.yaml_obj['up'][1], self.yaml_obj['up'][2],
           self.yaml_obj['width'], self.yaml_obj['height'], self.yaml_obj['field-of-view'],
           self.yaml_obj['focal-length'], self.yaml_obj['usteps'], self.yaml_obj['vsteps'])

class Light(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Light':
        if 'cache-size' not in obj:
            obj['cache-size'] = 65536
        if 'at' in obj:
            if 'to' in obj:
                if 'radius' in obj:
                    return CircleAreaLight.from_yaml(obj)
                return HemisphereLight.from_yaml(obj)
            return PointLight.from_yaml(obj)
        elif 'corner' in obj:
            return AreaLight.from_yaml(obj)

class HemisphereLight(Light):
    def __init__(self, yaml_obj):
        Light.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'HemisphereLight':
        return cls(obj)

    def c_repr(self, name):
        return """    /* hemisphere light {0} */
    Light hemisphere_light_{0} = all_lights + {0};
    Point hemisphere_light_{0}_point = {{ {1:.10f}, {2:.10f}, {3:.10f}, 1.0 }};
    Color hemisphere_light_{0}_intensity = color({4:.10f}, {5:.10f}, {6:.10f});
    Point hemisphere_light_{0}_to = {{ {7:.10f}, {8:.10f}, {9:.10f}, 0.0 }};
    hemisphere_light(hemisphere_light_{0}_point, hemisphere_light_{0}_to, hemisphere_light_{0}_intensity, hemisphere_light_{0});

    /* end point light {0} */
""".format(name,
           self.yaml_obj['at'][0], self.yaml_obj['at'][1], self.yaml_obj['at'][2],
           self.yaml_obj['intensity'][0], self.yaml_obj['intensity'][1], self.yaml_obj['intensity'][2],
           self.yaml_obj['to'][0], self.yaml_obj['to'][1], self.yaml_obj['to'][2])

class PointLight(Light):
    def __init__(self, yaml_obj):
        Light.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'PointLight':
        return cls(obj)

    def c_repr(self, name):
        return """    /* point light {0} */
    Light point_light_{0} = all_lights + {0};
    Point point_light_{0}_point = {{ {1:.10f}, {2:.10f}, {3:.10f}, 1.0 }};
    Color point_light_{0}_intensity = color({4:.10f}, {5:.10f}, {6:.10f});
    point_light(point_light_{0}_point, point_light_{0}_intensity, point_light_{0});

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
    Point area_light_{0}_corner = {{ {1:.10f}, {2:.10f}, {3:.10f}, 1.0}};
    Color area_light_{0}_intensity = color({4:.10f}, {5:.10f}, {6:.10f});
    Vector area_light_{0}_uvec = vector_init({7:.10f}, {8:.10f}, {9:.10f});
    Vector area_light_{0}_vvec = vector_init({10:.10f}, {11:.10f}, {12:.10f});
    area_light(area_light_{0}_corner, area_light_{0}_uvec, {13}/*usteps*/, area_light_{0}_vvec, {14}/*vsteps*/, {15}/*jitter*/, {16}/*cache_size*/, area_light_{0}_intensity, area_light_{0});

    /* end area light 0 */
""".format(name,
           self.yaml_obj['corner'][0], self.yaml_obj['corner'][1], self.yaml_obj['corner'][2],
           self.yaml_obj['intensity'][0], self.yaml_obj['intensity'][1], self.yaml_obj['intensity'][2],
           self.yaml_obj['uvec'][0], self.yaml_obj['uvec'][1], self.yaml_obj['uvec'][2],
           self.yaml_obj['vvec'][0], self.yaml_obj['vvec'][1], self.yaml_obj['vvec'][2],
           self.yaml_obj['usteps'], self.yaml_obj['vsteps'], bool_str, self.yaml_obj['cache-size'])


class CircleAreaLight(Light):
    def __init__(self, yaml_obj):
        Light.__init__(self, yaml_obj)

    @classmethod
    def from_yaml(cls, obj) -> 'CircleAreaLight':
        return cls(obj)

    def c_repr(self, name):
        bool_str = "false"
        if self.yaml_obj['jitter']:
            bool_str = "true"
        return """
    /* circle area light {0} */
    Light circle_area_light_{0} = all_lights + {0};
    Point circle_area_light_{0}_position = {{ {1:.10f}, {2:.10f}, {3:.10f}, 1.0}};
    Color circle_area_light_{0}_intensity = color({4:.10f}, {5:.10f}, {6:.10f});
    Point circle_area_light_{0}_to = vector_init({7:.10f}, {8:.10f}, {9:.10f});
    circle_light(circle_area_light_{0}_position, circle_area_light_{0}_to, {14:.10f}/*radius*/, {10}/*usteps*/, {11}/*vsteps*/, {12}/*jitter*/, {13}/*cache_size*/, circle_area_light_{0}_intensity, circle_area_light_{0});

    /* end circle area light 0 */
""".format(name,
           self.yaml_obj['at'][0], self.yaml_obj['at'][1], self.yaml_obj['at'][2],
           self.yaml_obj['intensity'][0], self.yaml_obj['intensity'][1], self.yaml_obj['intensity'][2],
           self.yaml_obj['to'][0], self.yaml_obj['to'][1], self.yaml_obj['to'][2],
           self.yaml_obj['usteps'], self.yaml_obj['vsteps'], bool_str, self.yaml_obj['cache-size'],
           self.yaml_obj['radius'])



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
