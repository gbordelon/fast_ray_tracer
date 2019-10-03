# fast_ray_tracer

Example usage:

```
python3 yaml_parser/yaml_parser.py scenes/bounding_boxes/bounding_boxes.yml >./main.c
make
date && time ./ray_tracer
convert /tmp/unclamped.ppm -depth 16 -quality 95 /tmp/unclamped.jpg
open /tmp/unclamped.ppm && open /tmp/unclamped.jpg
```
