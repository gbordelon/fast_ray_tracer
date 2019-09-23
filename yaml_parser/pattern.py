class Pattern(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Pattern':
        xform = matrix4x4identity()
        if 'transform' in obj:
            xform = Transform.from_yaml(obj['transform'])

        typ = obj['type']
        # base patterns
        if typ == 'checkers':
            color1 = color(*obj['colors'][0])
            color2 = color(*obj['colors'][1])
            return Checkers(color1=color1, color2=color2, transform=xform)
        elif typ == 'gradient':
            color1 = color(*obj['colors'][0])
            color2 = color(*obj['colors'][1])
            return Gradient(color1=color1, color2=color2, transform=xform)
        elif typ == 'radial-gradient':
            color1 = color(*obj['colors'][0])
            color2 = color(*obj['colors'][1])
            return RadialGradient(color1=color1, color2=color2, transform=xform)
        elif typ == 'rings':
            color1 = color(*obj['colors'][0])
            color2 = color(*obj['colors'][1])
            return Ring(color1=color1, color2=color2, transform=xform)
        elif typ == 'stripe':
            color1 = color(*obj['colors'][0])
            color2 = color(*obj['colors'][1])
            return Stripe(color1=color1, color2=color2, transform=xform)

        # nested patterns
        elif typ == 'blended':
            return BlendedPattern(pattern1=Pattern.from_yaml(obj['left']),
                                  pattern2=Pattern.from_yaml(obj['right']))
        elif typ == 'nested':
            return NestedPattern(pattern1=Pattern.from_yaml(obj['primary']),
                                 pattern2=Pattern.from_yaml(obj['left']),
                                 pattern3=Pattern.from_yaml(obj['right']))
        elif typ == 'perturbed':
            freq = 0.4
            scale_factor = 0.3
            octaves = 1
            if 'frequency' in obj:
                freq = np.float64(obj['frequency'])
            if 'scale-factor' in obj:
                scale_factor = np.float64(obj['scale-factor'])
            if 'octaves' in obj:
                octaves = obj['octaves']
            return PerturbedPattern(pattern1=Pattern.from_yaml(obj['primary']),
                                    frequency=freq,
                                    scale_factor=scale_factor,
                                    octaves=octaves)

        # uv mapped patterns
        elif typ == 'map':
            mapping = obj['mapping']
            if mapping == 'cube':
                return CubeMapPattern(uv_left=Pattern.uv_from_yaml(obj['left']),
                                      uv_front=Pattern.uv_from_yaml(obj['front']),
                                      uv_right=Pattern.uv_from_yaml(obj['right']),
                                      uv_back=Pattern.uv_from_yaml(obj['back']),
                                      uv_up=Pattern.uv_from_yaml(obj['up']),
                                      uv_down=Pattern.uv_from_yaml(obj['down']),
                                      uv_map=CubeUVMap)
            # TODO consider supporting both caps
            elif mapping == 'cylindrical':
                if 'uv_pattern' in obj:
                    uv_cap = Pattern.uv_from_yaml(obj['uv_pattern'])
                    uv_body = Pattern.uv_from_yaml(obj['uv_pattern'])
                else:
                    uv_cap = Pattern.uv_from_yaml(obj['top'])
                    uv_body = Pattern.uv_from_yaml(obj['front'])
                return CylinderMapPattern(uv_cap=uv_cap,
                                          uv_body=uv_body,
                                          uv_map=CylinderUVMap,
                                          transform=xform)
            elif mapping == 'planar':
                return TextureMapPattern(uv_pattern=Pattern.uv_from_yaml(obj['uv_pattern']),
                                         uv_map=PlaneUVMap,
                                         transform=xform)
            elif mapping == 'spherical':
                return TextureMapPattern(uv_pattern=Pattern.uv_from_yaml(obj['uv_pattern']),
                                         uv_map=SphereUVMap,
                                         transform=xform)
        raise ValueError('Unable to parse pattern type: {}'.format(typ))

    @classmethod
    def uv_from_yaml(cls, obj):
        typ = obj['type']
        if typ == 'checkers':
            # colors is a list
            color1 = color(*obj['colors'][0])
            color2 = color(*obj['colors'][1])
            return UVChecker(color1=color1, color2=color2, width=obj['width'], height=obj['height'])
        elif typ == 'align_check':
            colors = obj['colors'] # dict
            return UVAlignCheck(main_color=colors['main'],
                                upper_left=colors['ul'],
                                upper_right=colors['ur'],
                                bottom_left=colors['bl'],
                                bottom_right=colors['br'])
        elif typ == 'image':
            file_path = obj['file']
            return UVTexturePattern(tcanvas=TextureCanvas(file_path))

        raise ValueError('Unable to parse uv pattern type: {}'.format(typ))

    def c_repr(self, name):
        return "NOT IMPLEMENTED"
