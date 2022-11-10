# PPAR project 1: Earth elevation data #

`model_old.c` is the code that we were given unchanged
`model.c` solves the least square problem using the DGELS LAPACK function.

## Compiling ##

Just run the `Makefile` with the command `make`


## Calling a LAPACK function in C ##

Since LAPACK is written in fortran, some adjustments had to be made to call the DGELS function from C:
- At the beginning of model.c the we write the signature of the LAPACK function after the keyword `extern`. This is used by the compiler to link the function in the c file to the one in the library. The name of the function is the same as the original in fortran, but in lowercase.

    Example for the DGELS function: 
    ```C
    {
        extern void dgels(char * trans, int * m, int * n, int * nrhs, double * A, int * lda, double * B, int * ldb, double * work, int * lwork, int * info);
    }
    ```
- When calling the function, the name of the function should be written in lowercase (like in the signature) **followed by an undersore**.

    Example for the DGELS function:
    ```C
    {
        dgels_(&trans, &m, &n, &nrhs, A, &lda, data.V, &ldb, work, &lwork, &info, 0);
    }
    ```
    Since we're talking with a fortran function, **all arguments are passed by reference**.
    Depending on the fortran compiler, a `size_t` type argument may be required in the function. I think this indicates the size of the string argument of the function, if there is one. In this example there is none, so the argument is 0. This argument is passed by value.
- The `-llapack` and `-lblas` flags must be added to the `Makefile` to make the compiler link the library with the c file where the LAPACK functinos are being called
