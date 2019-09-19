CC = cc
CFLAGS = -std=c11 -pedantic -Wall
DEPS = linalg.h canvas.h perlin.h shapes.h renderer.h sphere.h plane.h cube.h cone.h cylinder.h
OBJ = linalg.o main.o canvas.o perlin.o shapes.o renderer.o sphere.o plane.o cube.o cone.o cylinder.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

ray_tracer: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

all: ray_tracer

clean:
	rm -f *.o ray_tracer
