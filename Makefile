LDLIBS = -llapacke -lm 
CFLAGS = -g -Wall -O3 

ALL: model_old validate

model_old: model_old.o harmonics.o
validate: validate.o harmonics.o 
model_old.o: harmonics.h
quality.o: harmonics.h
harmonics.o: harmonics.h

.PHONY: clean

clean:
	rm -f model_old validate *.o