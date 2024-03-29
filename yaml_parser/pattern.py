import subprocess
import os

from transform import Transform

class Pattern(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj) -> 'Pattern':
        return cls(yaml_obj=obj)

    def c_repr(self, name, resources):
        if len(self.yaml_obj) == 0:
            return "    Pattern pattern_{0} = NULL;\n".format(name)

        if 'transform' not in self.yaml_obj:
            self.yaml_obj['transform'] = []

        transform = Transform.from_yaml(self.yaml_obj['transform'])
        buf = "{}\n".format(transform.c_repr("pattern_{0}".format(name)))

        typ = self.yaml_obj['type']
        # base patterns
        if typ in ['checker', 'checkers', 'gradient', 'radial-gradient', 'rings', 'ring', 'stripe', 'stripes']:
            if typ == 'rings':
                typ = 'ring'
            if typ == 'checkers':
                typ = 'checker'
            if typ == 'stripes':
                typ = 'stripe'
            buf += """    Color pattern_{0}_color_0_raw = color({2:.10f}, {3:.10f}, {4:.10f});
    Color pattern_{0}_color_1_raw = color({5:.10f}, {6:.10f}, {7:.10f});
    Color pattern_{0}_color_0;
    Color pattern_{0}_color_1;
    color_space_fn(pattern_{0}_color_0_raw, pattern_{0}_color_0);
    color_space_fn(pattern_{0}_color_1_raw, pattern_{0}_color_1);
    Pattern pattern_{0} = {1}_pattern_alloc(pattern_{0}_color_0, pattern_{0}_color_1);

""".format(name, typ,
           self.yaml_obj['colors'][0][0], self.yaml_obj['colors'][0][1], self.yaml_obj['colors'][0][2],
           self.yaml_obj['colors'][1][0], self.yaml_obj['colors'][1][1], self.yaml_obj['colors'][1][2])

        # nested patterns
        elif typ == 'blended':
            pattern1 = Pattern.from_yaml(self.yaml_obj['left'])
            pattern2 = Pattern.from_yaml(self.yaml_obj['right'])

            buf += """    {1}
{2}
    Pattern pattern_{0} = blended_pattern_alloc(pattern_blended_{0}_0, pattern_blended_{0}_1);
""".format(name, pattern1.c_repr("blended_{0}_0".format(name), resources), pattern2.c_repr("blended_{0}_1".format(name), resources))
        elif typ == 'nested':
            pattern1 = Pattern.from_yaml(self.yaml_obj['primary'])
            pattern2 = Pattern.from_yaml(self.yaml_obj['left'])
            pattern3 = Pattern.from_yaml(self.yaml_obj['right'])

            buf += """    {1}
{2}
{3}
    Pattern pattern_{0} = nested_pattern_alloc(pattern_nested_{0}_0, pattern_nested_{0}_1, pattern_nested_{0}_2);
""".format(name, pattern1.c_repr("nested_{0}_0".format(name), resources), pattern2.c_repr("nested_{0}_1".format(name), resources), pattern1.c_repr("nested_{0}_2".format(name), resources))
        elif typ == 'perturbed':
            freq = 1.0
            scale_factor = 0.01
            persistence = 0.7
            octaves = 1
            seed = 0
            if 'frequency' in self.yaml_obj:
                freq = self.yaml_obj['frequency']
            if 'scale-factor' in self.yaml_obj:
                scale_factor = self.yaml_obj['scale-factor']
            if 'octaves' in self.yaml_obj:
                octaves = self.yaml_obj['octaves']
            if 'persistence' in self.yaml_obj:
                persistence = self.yaml_obj['persistence']
            if 'seed' in self.yaml_obj:
                seed = self.yaml_obj['seed']

            pattern1 = Pattern.from_yaml(self.yaml_obj['primary'])

            buf += """    {1}
    Pattern pattern_{0} = perturbed_pattern_alloc(pattern_perturbed_{0}_0, {2:.10f}, {3:.10f}, {4:.10f}, {5}, {6});
""".format(name, pattern1.c_repr("perturbed_{0}_0".format(name), resources), freq, scale_factor, persistence, octaves, seed)

        # uv mapped patterns
        elif typ == 'map':
            mapping = self.yaml_obj['mapping']
            if mapping in ['cube', 'cubic']:
                uv_left = UVPattern.uv_from_yaml(self.yaml_obj['left'])
                uv_front = UVPattern.uv_from_yaml(self.yaml_obj['front'])
                uv_right = UVPattern.uv_from_yaml(self.yaml_obj['right'])
                uv_back = UVPattern.uv_from_yaml(self.yaml_obj['back'])
                uv_up = UVPattern.uv_from_yaml(self.yaml_obj['up'])
                uv_down = UVPattern.uv_from_yaml(self.yaml_obj['down'])

                buf += """    Pattern pattern_{0} = array_of_patterns(7);
    Pattern pattern_{0}_right = pattern_{0} + 1;
    Pattern pattern_{0}_left = pattern_{0} + 2;
    Pattern pattern_{0}_up = pattern_{0} + 3;
    Pattern pattern_{0}_down = pattern_{0} + 4;
    Pattern pattern_{0}_front = pattern_{0} + 5;
    Pattern pattern_{0}_back = pattern_{0} + 6;

{1}
{2}
{3}
{4}
{5}
{6}

    texture_map_pattern(pattern_{0}_right, CUBE_UV_MAP, pattern_{0});
""".format(name,
           uv_right.c_repr('pattern_{0}_right'.format(name), resources),
           uv_left.c_repr('pattern_{0}_left'.format(name), resources),
           uv_up.c_repr('pattern_{0}_up'.format(name), resources),
           uv_down.c_repr('pattern_{0}_down'.format(name), resources),
           uv_front.c_repr('pattern_{0}_front'.format(name), resources),
           uv_back.c_repr('pattern_{0}_back'.format(name), resources))

            elif mapping in ['cylindrical', 'cylinder']:
                # handle single pattern and three patterns for cylinder
                if 'uv_pattern' in self.yaml_obj:
                    uv_top = UVPattern.uv_from_yaml(self.yaml_obj['uv_pattern'])
                    uv_bottom = uv_top
                    uv_body = uv_top
                else:
                    uv_top = UVPattern.uv_from_yaml(self.yaml_obj['top'])
                    uv_bottom = UVPattern.uv_from_yaml(self.yaml_obj['bottom'])
                    uv_body = UVPattern.uv_from_yaml(self.yaml_obj['front'])

                buf += """    Pattern pattern_{0} = array_of_patterns(4);
    Pattern pattern_{0}_body = pattern_{0} + 1;
    Pattern pattern_{0}_top = pattern_{0} + 2;
    Pattern pattern_{0}_bottom = pattern_{0} + 3;
""".format(name)

                buf += """{1}
{2}
{3}

    texture_map_pattern(pattern_{0}_body, CYLINDER_UV_MAP, pattern_{0});
""".format(name,
           uv_body.c_repr('pattern_{0}_body'.format(name), resources),
           uv_top.c_repr('pattern_{0}_top'.format(name), resources),
           uv_bottom.c_repr('pattern_{0}_bottom'.format(name), resources))

            elif mapping in ['triangular', 'triangle']:
                uv_pattern = UVPattern.uv_from_yaml(self.yaml_obj['uv_pattern'])
                buf += """    Pattern pattern_{0} = array_of_patterns(2);
    Pattern pattern_{0}_body = pattern_{0} + 1;
{1}
    texture_map_pattern(pattern_{0}_body, TRIANGLE_UV_MAP, pattern_{0});
""".format(name, uv_pattern.c_repr('pattern_{0}_body'.format(name), resources))
            elif mapping in ['planar', 'plane']:
                uv_pattern = UVPattern.uv_from_yaml(self.yaml_obj['uv_pattern'])
                buf += """    Pattern pattern_{0} = array_of_patterns(2);
    Pattern pattern_{0}_body = pattern_{0} + 1;
{1}
    texture_map_pattern(pattern_{0}_body, PLANE_UV_MAP, pattern_{0});
""".format(name, uv_pattern.c_repr('pattern_{0}_body'.format(name), resources))
            elif mapping in ['spherical', 'sphere']:
                uv_pattern = UVPattern.uv_from_yaml(self.yaml_obj['uv_pattern'])
                buf += """    Pattern pattern_{0} = array_of_patterns(2);
    Pattern pattern_{0}_body = pattern_{0} + 1;
{1}
    texture_map_pattern(pattern_{0}_body, SPHERE_UV_MAP, pattern_{0});
""".format(name, uv_pattern.c_repr('pattern_{0}_body'.format(name), resources))
            elif mapping in ['toroidal', 'toroid', 'torus']:
                uv_pattern = UVPattern.uv_from_yaml(self.yaml_obj['uv_pattern'])
                buf += """    Pattern pattern_{0} = array_of_patterns(2);
    Pattern pattern_{0}_body = pattern_{0} + 1;
{1}
    texture_map_pattern(pattern_{0}_body, TOROID_UV_MAP, pattern_{0});
""".format(name, uv_pattern.c_repr('pattern_{0}_body'.format(name), resources))
        else:
            raise ValueError('Unable to parse pattern type: {}'.format(typ))

        buf += "    pattern_set_transform(pattern_{0}, transform_pattern_{0});\n".format(name)
        return buf


class UVPattern(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def uv_from_yaml(cls, obj):
        return cls(obj)

    def c_repr(self, name, resources):
        typ = self.yaml_obj['type']
        if typ in ['checkers', 'check']:
            buf = """    Color {0}_color_0_raw = color({2:.10f}, {3:.10f}, {4:.10f});
    Color {0}_color_1_raw = color({5:.10f}, {6:.10f}, {7:.10f});
    Color {0}_color_0;
    Color {0}_color_1;
    color_space_fn({0}_color_0_raw, {0}_color_0);
    color_space_fn({0}_color_1_raw, {0}_color_1);
    uv_check_pattern({0}_color_0, {0}_color_1, {8}, {9}, {0});

""".format(name,
           "",
           self.yaml_obj['colors'][0][0],
           self.yaml_obj['colors'][0][1],
           self.yaml_obj['colors'][0][2],
           self.yaml_obj['colors'][1][0],
           self.yaml_obj['colors'][1][1],
           self.yaml_obj['colors'][1][2],
           self.yaml_obj['width'],
           self.yaml_obj['height'])

        elif typ in ['align_check', 'align-check']:
            if typ == 'align-check':
                typ = 'align_check'
            colors = self.yaml_obj['colors'] # dict

            buf = """    Color {0}_color_0_raw = color({2:.10f}, {3:.10f}, {4:.10f});
    Color {0}_color_1_raw = color({5:.10f}, {6:.10f}, {7:.10f});
    Color {0}_color_2_raw = color({8:.10f}, {9:.10f}, {10:.10f});
    Color {0}_color_3_raw = color({11:.10f}, {12:.10f}, {13:.10f});
    Color {0}_color_4_raw = color({14:.10f}, {15:.10f}, {16:.10f});
    Color {0}_color_0;
    Color {0}_color_1;
    Color {0}_color_2;
    Color {0}_color_3;
    Color {0}_color_4;
    color_space_fn({0}_color_0_raw, {0}_color_0);
    color_space_fn({0}_color_1_raw, {0}_color_1);
    color_space_fn({0}_color_2_raw, {0}_color_2);
    color_space_fn({0}_color_3_raw, {0}_color_3);
    color_space_fn({0}_color_4_raw, {0}_color_4);
    uv_align_check_pattern({0}_color_0, {0}_color_1, {0}_color_2, {0}_color_3, {0}_color_4, {0});

""".format(name,
           "",
           self.yaml_obj['colors']['main'][0],
           self.yaml_obj['colors']['main'][1],
           self.yaml_obj['colors']['main'][2],
           self.yaml_obj['colors']['ul'][0],
           self.yaml_obj['colors']['ul'][1],
           self.yaml_obj['colors']['ul'][2],
           self.yaml_obj['colors']['ur'][0],
           self.yaml_obj['colors']['ur'][1],
           self.yaml_obj['colors']['ur'][2],
           self.yaml_obj['colors']['bl'][0],
           self.yaml_obj['colors']['bl'][1],
           self.yaml_obj['colors']['bl'][2],
           self.yaml_obj['colors']['br'][0],
           self.yaml_obj['colors']['br'][1],
           self.yaml_obj['colors']['br'][2])
        elif typ == 'image':
            file_path = self.yaml_obj['file']
            png_file_path = file_path[:-3] + 'png';
            if file_path[-3:] != 'png':
                # check for existence of ppm extension file name
                if not (os.path.exists(png_file_path) and os.path.isfile(png_file_path)):
                    # if it does not exist, create it with 'convert'
                    subprocess.run(['convert', file_path, '-compress', 'none', '-quality', '95', png_file_path])

            if file_path not in resources:
                uv_pattern_name = 'pattern_{}'.format(file_path.replace('/', '_').replace('-','_').replace('.','_'))
                color_space_fn_name = 'rgb_to_rgb'
                if name.find('Ka') > 0 or name.find('Kd') > 0:
                    color_space_fn_name = 'color_space_fn'
                buf = """    if (access("{1}", F_OK ) == -1 ) {{
        printf("file '{1}' does not exist.");
        return 1;
    }}
    printf("Loading resource '{1}'... ");
    fflush(stdout);
    Canvas {2};
    read_png(&{2}, "{1}", false, {3});
    uv_texture_pattern({2}, {0});
    printf("Done!\\n");
    fflush(stdout);
""".format(name, png_file_path, uv_pattern_name, color_space_fn_name)
                resources[file_path] = uv_pattern_name
            else:
                buf = """    uv_texture_pattern({1}, {0});
""".format(name, resources[file_path])
        else:
            raise ValueError('Unable to parse uv pattern type: {}'.format(typ))

        return buf
