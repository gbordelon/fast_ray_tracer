CC = cc
CFLAGS = -std=c99 -pedantic -Wall
DEPS = linalg.h canvas.h perlin.h
OBJ = linalg.o main.o canvas.o perlin.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

ray_tracer: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

all: ray_tracer

clean:
	rm -f *.o ray_tracer
