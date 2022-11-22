CC = mpicc
LDLIBS = -lm -pg
CFLAGS = -g -Wall -O3 

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi)
LDLIBS += $(shell pkg-config --libs scalapack-openmpi)

ALL: pmodel2 validate

pmodel2: pmodel2.o harmonics.o	
validate: validate.o harmonics.o 
pmodel2.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f pmodel2 validate *.o