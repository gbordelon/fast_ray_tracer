from shapes import Group

class OBJParser(object):
    @classmethod
    def from_yaml(cls, obj, defines) -> 'Group':
        return cls(obj, defines)

    def __init__(self, yaml_obj, defines=None):
        self.yaml_obj = yaml_obj
        self.defines = defines
        self.vertices = []
        self.normals = []
        self.textures = []
        self.vertices.append(None)
        self.normals.append(None)
        self.named_groups = {'##default_group' : {'add': 'group', 'children':[]}}
        self.current_group_name = '##default_group'
        self.default_group = self.named_groups['##default_group']
        self._parse_file()

    def _parse_file(self):
        with open(self.yaml_obj['file'], 'r') as f:
            for line in f:
                _line = line.strip()
                if len(_line) > 0:
                    self._parse_line(_line)

    def _parse_line(self, line):
        lsp = line.split(' ')
        if '' in lsp:
            lsp.remove('')
        if lsp[0] == 'v':
            self.vertices.append(self._parse_vertex(lsp))
        elif lsp[0] == 'vn':
            self.normals.append(self._parse_normal(lsp))
        elif lsp[0] == 'vt':
            self.textures.append(self._parse_texture(lsp))
        elif lsp[0] == 'f':
            self.named_groups[self.current_group_name]['children'].extend(self._fan_triangulation(lsp))
        elif lsp[0] == 'g':
            self._parse_group(lsp)

    def _parse_group(self, line):
        if line[1] not in self.named_groups:
            self.named_groups[line[1]] = {'add':'group', 'children':[]}
        self.current_group_name = line[1]

    def _parse_normal(self, line):
        try:
            return [float(line[1]), float(line[2]), float(line[3])]
        except ValueError as e:
            print(e, line)
            raise

    def _parse_vertex(self, line):
        try:
            return [float(line[1]), float(line[2]), float(line[3])]
        except ValueError as e:
            print(e, line)
            raise

    def _parse_texture(self, line):
        try:
            return [float(line[1]), float(line[2]), float(line[3])]
        except ValueError as e:
            print(e, line)
            raise

    # TODO should also check for existence of // instead of just /
    def _fan_triangulation(self, line):
        triangles = []
        if '/' in line[1]:
            v1 = line[1].split('/')[0]
            t1 = line[1].split('/')[1]
            n1 = line[1].split('/')[2]
        else:
            v1 = line[1]

        for i in range(2, len(line) - 1):
            if '/' in line[i]:
                v2 = line[i].split('/')[0]
                t2 = line[i].split('/')[1]
                n2 = line[i].split('/')[2]
            else:
                v2 = line[i]

            if '/' in line[i+1]:
                v3 = line[i+1].split('/')[0]
                t3 = line[i+1].split('/')[1]
                n3 = line[i+1].split('/')[2]
            else:
                v3 = line[i+1]

            try:
                if len(self.normals) > 0:
                    tri =  {'add':'smooth-triangle',
                            'p1':self.vertices[int(v1)],
                            'p2':self.vertices[int(v2)],
                            'p3':self.vertices[int(v3)],
                            'n1':self.normals[int(n1)],
                            'n2':self.normals[int(n2)],
                            'n3':self.normals[int(n3)]}
                else:
                    tri =  {'add':'triangle',
                            'p1':self.vertices[int(v1)],
                            'p2':self.vertices[int(v2)],
                            'p3':self.vertices[int(v3)]}

                if len(self.textures) > 0:
                    tri['t1'] = self.textures[int(t1)]
                    tri['t2'] = self.textures[int(t2)]
                    tri['t3'] = self.textures[int(t3)]

                triangles.append(tri)
            except IndexError as e:
                print(e, v1, v2, v3, n1, n2, n3, t1, t2, t3)
                raise

        return triangles

    def c_repr(self, name, parent_name, offset):
        obj = self.yaml_obj
        if 'material' not in obj:
            obj['material'] = {}

        if 'transform' not in obj:
            obj['transform'] = []

        top_group = {'add':'group',
                     'children': [v for k,v in self.named_groups.items() if len(v['children']) > 0]}
        top_group['material'] = obj['material']
        top_group['transform'] = obj['transform']

        item = Group.from_yaml(top_group, self.defines)
        return item.c_repr(name, parent_name, offset)
        

def obj_to_group(parser):
    return obj_group

