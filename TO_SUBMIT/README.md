# PPAR project 1: Earth elevation data #
by Marco PARADINA, Diego ARBULU SEDANO, Vignesh ANANTHARAMAKRISHNAN
<br/><br/>

`pmodel.c` contains the ScaLAPACK implementation and computes the least square solution

`scalapack.h` is the header file for all the ScaLAPACK and Cblacs routines' signatures

`install_scalapack` is a bash script to install ScaLAPACK on all the nodes that are being used for a job on grid5000

## Installing ScaLAPACK ##

Before compiling, the ScaLAPACK library must be installed. On grid5000 this can be done using the `install_scalapack.sh` bash script. From one of the nodes that are being used for the computation, use this command:

`bash install_scalapack.sh`

## Setting parameters ##

* The dimension of the process grid must be set accordingly to how many processes are being used for the computation. These can be set in the `NPROW`, `NPCOL` macros in `pmodel.c`.
* The blocking parameters can be modified. This can be done at lines 128-131 of `pmodel.c` by changing the values of `mb_A`, `nb_A`, `mb_b`, `nb_b` variables.


## Compilation ##

Run the `Makefile` with the command `make`

## Execution on grid5000 ##

To launch the computation on grid5000 use this command:

`mpirun -n [number of processes] -machinefile $OAR_NODEFILE ./pmodel --data [input file name] --model [output file name] --npoint [number of points] --lmax [lmax]` 

## Validation ##

To validate the result obtained from the computation, use `validate.c` as follows:

`./validate --data [input file name] --model [output file name] --lmax [lmax used in the computation]`



`pmodel.c` produces 2 outputs:
* The file containing the least-square solution
* A file containig the residual sum of squares after solving the least-square-problem and the time it took to execute, called `timer_file.txt`
