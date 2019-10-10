# fast_ray_tracer

Example usage:

```
python3 yaml_parser/yaml_parser.py scenes/bounding_boxes/bounding_boxes.yml >./main.c
make
date
time ./ray_tracer
open /tmp/srgb_render.ppm
open /tmp/srgb_render.png
```

TODO
* Refactor yaml parser logic to not load a file multiple times but use shape_copy on the parent group for object files.
* Refactor canvas/texture maps to only keep one image in memory even if multiple patterns try to load the file.
    * obj_loader default vertex num, default group num
    * ~~light source cache size~~
    * ~~Add the ability to configure color space in a yaml config file~~
* Refactor patterns into multiple source files
* Refactor materials out of shapes.{c,h}
* Add a circular area light
* Add a spot light
* Use realloc instead of malloc when resizing arrays for intersections and group.children
* Implement triangle mesh support
* Refactor photon mapping code in light.c to use CMJ instead of random sampling for photon emission
* Finish adding HSL color support
* Remove malloc/free from the photon map where possible.
* Implement projection map for generating the caustics photon map.
* Implement selective projection maps. Perhaps object opt-in.
* Genericize patterns so they don't use the Color type
* Add const keyword to function signatures where appropriate.
* Investigate BLAS/LAPACK
* Add heirarchical yaml parsing so I can have a global config which is read first before a scene's yaml file.
* Refactor shapes such that one only needs to include shapes.h
* Implement MTL file parsing.
* Add parameters to the config parser
    * epsilon
* ~~Implement photon tracing~~
* ~~Implement photon mapping for caustics~~
    * Balance radiance estimate
* ~~Implement photon mapping for diffuse reflection~~
    * Balance radiance estimate
* ~~Improve correlated multi-jitter logic~~
* ~~Refactor sampling and jittering from camera, renderer, lights into src/libs/sampler/sampler.{c,h}~~
* ~~Refactor core-selection code from renderer.c into a library~~
* ~~Support a thread pool so rows can be handed out independently.~~
* ~~Add GI parameters to the config parser.~~
    * ~~photon map cone constant~~
