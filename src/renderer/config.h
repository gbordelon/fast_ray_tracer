#ifndef CONFIG_H
#define CONFIG_H

struct direct_illumination {
    bool include_ambient;
    bool include_diffuse;
    bool include_specular_highlight;
    bool include_specular;
};

struct global_illumination {
    bool include_caustics;
    bool include_final_gather;
    size_t usteps;
    size_t vsteps;
    size_t irradiance_estimate_num;
    double irradiance_estimate_radius;
    size_t photon_count;
};

struct illumination_config {
    bool include_direct;
    bool include_global;
    struct direct_illumination di;
    struct global_illumination gi;
    // These two are mutually exclusive and will override other settings
    bool debug_visualize_photon_map;
    bool debug_visualize_soft_indirect;
};

struct threading_config {
    size_t num_threads;
};

struct scene_config {
    size_t divide_threshold;
};

typedef struct global_config {
    struct illumination_config illumination;
    struct threading_config threading;
    struct scene_config scene;

} *Global_config;

#endif
