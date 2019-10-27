#ifndef CONFIG_H
#define CONFIG_H

struct direct_illumination {
    bool include_ambient;
    bool include_diffuse;
    bool include_specular_highlight;
    bool include_specular;
    size_t path_length;
};

struct global_illumination {
    bool include_caustics;
    bool include_final_gather;
    size_t usteps;
    size_t vsteps;
    size_t irradiance_estimate_num;
    double irradiance_estimate_radius;
    double irradiance_estimate_cone_filter_k;
    size_t photon_count;
    size_t path_length;
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

enum color_space_type {
    SRGB,
    RGB,
    HSL,
    XYZ,
    XYY,
    LAB
};

struct output_config {
    const char *file_path;
    enum color_space_type color_space;
};

typedef struct global_config {
    struct illumination_config illumination;
    struct threading_config threading;
    struct scene_config scene;
    struct output_config output;

} *Global_config;

#endif
