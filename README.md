# PPAR project #1

secret bill gates project to destroy the world by computing the elevation of places

model_old.c is the code that we were given unchanged
model.c solves the least square problem using the DGELS LAPACK function.

Since LAPACK is written in fortran, some adjustments had to be made to call the DGELS function from C:
-at the beginning of model.c the signature for DGELS is states. This is used by the compiler to link the function in the c file to the one in the library.
-Since we're talking with a fortran functionl, all parameters are passed by reference to dgels_
-The -llapack flag was added to the Makefile to make the compiler link the library with model.c