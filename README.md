# fast_ray_tracer

Example usage:

```
python3 yaml_parser/yaml_parser.py scenes/bounding_boxes/bounding_boxes.yml >./main.c
make
date && time ./ray_tracer
convert /tmp/unclamped.ppm -depth 16 -quality 95 /tmp/unclamped.jpg
open /tmp/unclamped.ppm && open /tmp/unclamped.jpg
```

TODO
* Implement triangle mesh support
* ~~Implement photon tracing~~
* Implement photon mapping for caustics
* Implement photon mapping for diffuse reflection
* ~~Improve correlated multi-jitter logic~~
* Refactor patterns into multiple source files
* Refactor shapes such that one only needs to include shapes.h
* Refactor materials out of shapes.{c,h}
* ~~Refactor sampling and jittering from camera, renderer, lights into src/libs/sampler/sampler.{c,h}~~
* Refactor core-selection code from renderer.c into a library
* Refactor yaml parser logic to not load a file multiple times but use shape_copy on the parent group for object files.
* Refactor canvas/texture maps to only keep one image in memory even if multiple patterns try to load the file.
* Investigate BLAS/LAPACK
* Add const keyword to function signatures where appropriate.
* Finish adding HSL color support
* Add the ability to configure color space in a yaml config file
* Add heirarchical yaml parsing so I can have a global config which is read first before a scene's yaml file.
* Add a spot light
* Genericize patterns so they don't use the Color type
* ~~Support a thread pool so rows can be handed out independently.~~
* Use realloc instead of malloc when resizing arrays for intersections and group.children
* Remove malloc/free from the photon map where possible.
* Implement MTL file parsing.
* Implement project map for generating the caustics photon map.
* Implement selective projection maps. Perhaps object opt-in.
* Add GI parameters to the config parser.
