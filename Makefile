CC = cc
CFLAGS = -std=c11 -pedantic -Wall -g
DEPS = linalg.h canvas.h perlin.h shapes.h renderer.h sphere.h plane.h cube.h cone.h cylinder.h triangle.h csg.h group.h bounding_box.h toroid.h Roots3and4.h obj_loader.h thpool.h
OBJ = linalg.o main.o canvas.o perlin.o shapes.o renderer.o sphere.o plane.o cube.o cone.o cylinder.o triangle.o csg.o group.o bounding_box.o toroid.o Roots3and4.c obj_loader.o thpool.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

ray_tracer: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

all: ray_tracer

clean:
	rm -f *.o ray_tracer
