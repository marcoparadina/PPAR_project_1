CC = mpicc
LDLIBS = -lm 
CFLAGS = -g -Wall -O3 

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi)	
LDLIBS += $(shell pkg-config --libs scalapack-openmpi)

ALL: pmodel validate

pmodel: pmodel.o harmonics.o	
validate: validate.o harmonics.o 
pmodel.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f pmodel validate *.o