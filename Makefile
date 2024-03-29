CC = mpicc
LDLIBS = -lm 
CFLAGS = -g -Wall -O3 

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi)	
LDLIBS += $(shell pkg-config --libs scalapack-openmpi)

ALL: pmodel3_clean validate

pmodel3_clean: pmodel3_clean.o harmonics.o	
validate: validate.o harmonics.o 
pmodel3_clean.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f pmodel3_clean validate *.o