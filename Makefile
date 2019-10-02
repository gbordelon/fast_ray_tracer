TOPDIR	:= .
SRCDIR	:= $(TOPDIR)/src
OBJDIR	:= $(TOPDIR)/obj
#VPATH	:= $(SRCDIR)/color:$(SRCDIR)/renderer:$(SRCDIR)/shapes:$(SRCDIR)/libs/canvas:$(SRCDIR)/libs/linalg:$(SRCDIR)/libs/obj_loader:$(SRCDIR)/libs/perlin:$(SRCDIR)/libs/photon_map:$(SRCDIR)/libs/quartic:$(SRCDIR)/libs/thpool
#VPATH	:= $(SRCDIR)

CFILES	:= c
OFILES	:= o
HFILES	:= h
CC		:= cc
CFLAGS	:= -std=c11 -pedantic -Wall -g
#-I$(SRCDIR)/color -I$(SRCDIR)/renderer -I$(SRCDIR)/shapes -I$(SRCDIR)/libs/canvas -I$(SRCDIR)/libs/linalg -I$(SRCDIR)/libs/obj_loader -I$(SRCDIR)/libs/perlin -I$(SRCDIR)/libs/photon_map -I$(SRCDIR)/libs/quartic:$(SRCDIR)/libs/thpool

C_NAMES	:= $(shell find . -name '*.c')
H_NAMES	:= $(shell find . -name '*.h')

EXE		:= $(TOPDIR)/ray_tracer
SOURCES	:= $(C_NAMES)
DEPS	:= $(H_NAMES)
OBJECTS	:= $(shell sed -e 's/.\/src/obj/g' <<< " $(C_NAMES:%.c=%.o) " )

MAIN	:= main.c

.PHONY: all alldefault clean exe

exe: $(EXE)


$(EXE): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJDIR)/%.$(OFILES): $(SRCDIR)/%.$(CFILES) $(DEPS)
	mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c -o $@ $<

$(SRCDIR)/%.$(CFILES):
	: $@
	: $<

clean:
	rm -f $(EXE) main.o
	rm -rf $(OBJDIR)

test:
	: $(OBJECTS)
