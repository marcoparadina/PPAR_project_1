CC = mpicc
LDLIBS = -lm -pg
CFLAGS = -g -Wall -O3 

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi)
LDLIBS += $(shell pkg-config --libs scalapack-openmpi)

ALL: model_old validate

model_old: model_old.o harmonics.o
validate: validate.o harmonics.o 
model_old.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f model_old validate *.o