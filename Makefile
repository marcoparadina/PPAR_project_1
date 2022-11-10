LDLIBS = -lm -llapack -lblas
CFLAGS = -g -Wall -O3 

ALL: model validate

model: model.o harmonics.o
validate: validate.o harmonics.o 
model.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f model validate *.o