LIB = libflexure.a

SRC_DIR = src

OBJECTS_LIB = \
flexure2D.o \
four1.o \
realft.o \
sinft.o

FLAGS = -Wall -c -O3 #-fbounds-check

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
