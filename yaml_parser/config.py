class GlobalConfig(object):
    def __init__(self, yaml_obj):
        self.yaml_obj = yaml_obj

    @classmethod
    def from_yaml(cls, obj):
        if 'thread-count' not in obj:
            obj['thread-count'] = 4
        if 'timeout' not in obj:
            obj['timeout'] = 30
        if 'divide-threshold' not in obj:
            obj['divide-threshold'] = 1
        if 'clamping' not in obj:
            obj['clamping'] = False
        if 'photon-count' not in obj:
            obj['photon-count'] = 10000

        return cls(yaml_obj=obj)

    def c_repr(self):
        bool_str = "false"
        if self.yaml_obj['clamping']:
            bool_str = "true"
        return """    /* config */
    size_t thread_count = {0};
//    size_t timeout = {1};
    size_t divide_threshold = {2};
//    bool clamping = {3};
    size_t photon_count = {4};
    /* end config */
""".format(self.yaml_obj['thread-count'],
           self.yaml_obj['timeout'],
           self.yaml_obj['divide-threshold'],
           bool_str,
           self.yaml_obj['photon-count'])
