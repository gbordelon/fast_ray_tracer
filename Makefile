CC = cc
CFLAGS = -std=c99 -pedantic -Wall
DEPS = linalg.h canvas.h
OBJ = linalg.o main.o canvas.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

ray_tracer: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

all: ray_tracer

clean:
	rm -f *.o ray_tracer
