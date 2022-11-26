CC = mpicc
LDLIBS = -lm -pg
CFLAGS = -g -Wall -O3 

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi)	
LDLIBS += $(shell pkg-config --libs scalapack-openmpi)

ALL: pmodel3 validate

pmodel3: pmodel3.o harmonics.o	
validate: validate.o harmonics.o 
pmodel3.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f pmodel3 validate *.o