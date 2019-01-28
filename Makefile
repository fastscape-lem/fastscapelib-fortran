LIB = libfastscape.a

SRC_DIR = src

OBJECTS_LIB = \
FastScape_ctx.o \
Advect.o \
Diffusion.o \
LocalMinima.o \
StreamPowerLaw.o \
Strati.o \
VTK.o \
FastScape_api.o

FLAGS = -c -O3 #-fbounds-check

.PHONY: all clean

all: $(LIB)

$(LIB): $(OBJECTS_LIB)
	libtool -static -o $@ $^

%.o: $(SRC_DIR)/%.f90
	gfortran $(FLAGS) -c $< -o $@

%.o: $(SRC_DIR)/%.f
	gfortran $(FLAGS) -c $< -o $@

%.o: $(SRC_DIR)/%.c
	gcc $(FLAGS)-c $< -o $@

clean:
	$(RM) *.o
