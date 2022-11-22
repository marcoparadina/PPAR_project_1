#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <lapack.h>
#include "harmonics.h"
#include "mpi.h"
#include "scalapack.h"	//not sure this is the right library to import, might have to change

/* Signature of SCALAPACK function */
extern void pdgels(char * trans, int * m, int * n, int * nrhs, double * A, int * lda, double * B, int * ldb, double * work, int * lwork, int * info);

int lmax = -1;
int npoint;
char * data_filename;
char * model_filename;

void usage(char ** argv)
{
	printf("%s [OPTIONS]\n\n", argv[0]);
	printf("Options:\n");
	printf("--data FILENAME              input file containing experimental data points\n");
	printf("--model FILENAME             output file containing the model\n");
	printf("--npoint N                   number of points to read\n");
	printf("--lmax N                     order of the model\n");
	printf("\n");
	exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
	struct option longopts[5] = {
		{"data", required_argument, NULL, 'd'},
		{"npoint", required_argument, NULL, 'n'},
		{"lmax", required_argument, NULL, 'l'},
		{"model", required_argument, NULL, 'm'},
		{NULL, 0, NULL, 0}
	};
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'd':
			data_filename = optarg;
			break;
		case 'm':
			model_filename = optarg;
			break;
		case 'n':
			npoint = atoi(optarg);
			break;
		case 'l':
			lmax = atoll(optarg);
			break;
		default:
			errx(1, "Unknown option\n");
		}
	}
	/* missing required args? */
	if (data_filename == NULL || model_filename == NULL || lmax < 0 || npoint <= 0)
		usage(argv);
}

int main(int argc, char ** argv){

    MPI_Init(NULL, NULL);
    int sysctx; 

    // MPI parameters 
	int rank = 0; //The rank of this process
	int p = 0; //Number of processes

    //Set up BLACS here. How? I don't fucking know

    if(rank == 0){
        process_command_line_options(argc, argv);

        /* preparations and memory allocation */
        int nvar = (lmax + 1) * (lmax + 1);
        printf("Linear Least Squares with dimension %d x %d\n", npoint, nvar);
        if (nvar > npoint)
            errx(1, "not enough data points");

        long matrix_size = sizeof(double) * nvar * npoint;
        char hsize[16];
        human_format(hsize, matrix_size);
        printf("Matrix size: %sB\n", hsize);

        double *A = malloc(matrix_size);
        if (A == NULL)
            err(1, "cannot allocate matrix");

        double * P = malloc((lmax + 1) * (lmax + 1) * sizeof(*P));
        double * v = malloc(npoint * sizeof(*v));
        if (P == NULL || v == NULL)
            err(1, "cannot allocate data points\n");

        printf("Reading data points from %s\n", data_filename);
        struct data_points data;
        load_data_points(data_filename, npoint, &data);
        printf("Successfully read %d data points\n", npoint);
        
        printf("Building matrix\n");
        struct spherical_harmonics model;
        setup_spherical_harmonics(lmax, &model);

        for (int i = 0; i < npoint; i++) {
            computeP(&model, P, sin(data.phi[i]));
            
            for (int l = 0; l <= lmax; l++) {
                /* zonal term */
                A[i + npoint * CT(l, 0)] = P[PT(l, 0)];
        
                /* tesseral terms */
                for (int m = 1; m <= l; m++) {
                    A[i + npoint * CT(l, m)] = P[PT(l, m)] * cos(m * data.lambda[i]);
                    A[i + npoint * ST(l, m)] = P[PT(l, m)] * sin(m * data.lambda[i]);
                }
            }
        }
        
        double FLOP = 2. * nvar * nvar * npoint;
        char hflop[16];
        human_format(hflop, FLOP);
        printf("Least Squares (%sFLOP)\n", hflop);
        double start = wtime();
    }

	else{
        //call pdgels on the submatrices of this process
        int m, n, mb, nb, nrhs, lda, ldb, lwork, info;
        m = npoint;	//Number of rows of the A matrix
        n = nvar;	//Number of columns of the A matrix
        mb = 1;  //Number if rows of a 2D block
        nb = 1; //Number if columns of a 2D block
        nrhs = 1;	//Number of columns of the x and B matrices, where ||Ax-B|| is the system that we're minimizing
        lda = npoint; //Leading dimension of the A matrix, i.e. its number of columns
        ldb = npoint; //Leading dimension of the B matrix, i.e. its number of columns
        lwork = 4*nvar;
        info = 0; 
        double work[4*nvar];
        int descA[9], descB[9];
	
    //Initialize descriptors for global matrices
    descinit_(descA, &ctx, &npoint, &nvar, &mb, &nb, &i_zero, &i_zero, )

    }

	
    /* Define the parameters to give to pdgels */
    

	/*************************LAPACK CALL***************************/

	
	lapack_int m, n, nrhs, lda, ldb, lwork, info;
	char trans = 'N';
	m = npoint;	//Number of rows of the A matrix
	n = nvar;	//Number of columns of the A matrix
	nrhs = 1;	//Number of columns of the x and B matrices, where ||Ax-B|| is the system that we're minimizing
	lda = npoint; //Leading dimension of the A matrix, i.e. its number of columns
	ldb = npoint; //Leading dimension of the B matrix, i.e. its number of columns
	lwork = 4*nvar;
	info = 0; 
	double work[4*nvar];

	/*Notice how here dgels is followed by and underscore. Why? No idea, but it doesn't work otherwise.
	  The last parameter is a size_t. I think it represents the size of the string argument, if there is one. Here there's no such 
	  argument, so I set it to 0, but any value seems to work. This is used by the fortran compiler in some way, and NEEDS to be declared,
	  otherwise an error message "too few arguments" will be given. 
	*/

	dgels_(&trans, &m, &n, &nrhs, A, &lda, data.V, &ldb, work, &lwork, &info, 0);
	
	
	/**************************************************************/
	
	double t = wtime()  - start;
	double FLOPS = FLOP / t;
	char hflops[16];
	human_format(hflops, FLOPS);
	printf("Completed in %.1f s (%s FLOPS)\n", t, hflops);

	double res = 0;
	for (int j = nvar; j < npoint; j++)
		res += data.V[j] * data.V[j];
	printf("residual sum of squares %g\n", res);

	printf("Saving model in %s\n", model_filename);
	FILE *g = fopen(model_filename, "w");
	if (g == NULL)
		err(1, "cannot open %s for writing\n", model_filename);
	for (int l = 0; l <= lmax; l++) {
		fprintf(g, "%d\t0\t%.18g\t0\n", l, data.V[CT(l, 0)]);
		for (int m = 1; m <= l; m++)
			fprintf(g, "%d\t%d\t%.18g\t%.18g\n", l, m, data.V[CT(l, m)], data.V[ST(l, m)]);
	}
	return(info);
}

int main(int argc, char** argv){    

    // Useful constants for ScaLapack
    const MKL_INT i_one = 1, i_negone = -1, i_zero = 0;
    double zero=0.0E+0, one=1.0E+0;

    //Used for BLACS grid
    int nprow=2, npcol=2, myrow=-1, mycol=-1;

    //dimension for global matrix, Row and col blocking params
    int M, N, mb, nb, LWORK, INFO=0, ictxt=0;
    int descA_distr[9],descA[9];
    int lld=0;
    double *A;

    //dimension for local matrix
    int mp, nq, lld_distr;
    double *A_distr;
    double *TAU,*WORK;

    //counters, seeds
    int i=0, itemp, seed;

    //call MPI init commands
    MPI_Init(NULL, NULL);

    // Input your parameters: m, n - matrix dimensions, mb, nb - blocking parameters,
    // nprow, npcol - grid dimensions
    //Want 2x2 grid, will change in future
    nprow = 2;
    npcol = 2;
    //Matrix Size
    M=16;
    N=16;
    //Matrix blocks
    mb=2;
    nb=2;
    //for local work
    LWORK=nb*M;

    // Part with invoking of ScaLAPACK routines. Initialize process grid, first
    /*
    //Cblacs routines, do not work
    Cblacs_pinfo( &iam, &nprocs ) ;
    Cblacs_get( -1, 0, &ictxt );
    Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

    printf("Thread=%d, Cblacs_pinfo, iam=%d,nprocs=%d\n",myrank_mpi,iam,nprocs);
    printf("Thread=%d, After info, nprow=%d, npcol=%d, myrow=%d, mycol=%d\n",myrank_mpi,nprow,npcol,myrow,mycol);
    */

    //blacs Method, does not work
    printf("Thread %d: Before init, nprow=%d, npcol=%d, myrow=%d, mycol=%d\n",myrank_mpi,nprow,npcol,myrow,mycol);
    //printf("before, ictxt=%d\n",ictxt);
    blacs_get_( &i_zero, &i_zero, &ictxt );
    blacs_gridinit_( &ictxt, "R", &nprow, &npcol );
    blacs_gridinfo_( &ictxt, &nprow, &npcol, &myrow, &mycol );
    printf("Thread %d: After init, nprow=%d, npcol=%d, myrow=%d, mycol=%d\n",myrank_mpi,nprow,npcol,myrow,mycol);
    //Generate Data
    if ( myrow==0 && mycol==0 ){
    A = malloc(M*N*sizeof(double));
    //just generate some random samples for now,
    for(i = 0; i <M*N; i++) {
    A[i]=rand()%1000000;
    }
    //FOR DEBUG
    /*A[0]=1.0;
    A[1]=2.0;
    A[2]=3.0;
    A[3]=4.0;
    A[4]=4.0;
    A[5]=3.0;
    A[6]=2.0;
    A[7]=1.0;
    A[8]=1.0;
    A[9]=2.0;
    A[10]=3.0;
    A[11]=4.0;
    A[12]=4.0;
    A[13]=3.0;
    A[14]=2.0;
    A[15]=1.0;*/
    }else{
    A = NULL;
    //other processes don't contain parts of A
    }

    // Compute dimensions of local part of distributed matrix A_distr
    /*
    * DEBUG
    printf("M=%d\n",M);
    printf("mb=%d\n",mb);
    printf("myrow=%d\n",myrow);
    printf("izero=%d\n",i_zero);
    printf("nprow=%d\n",nprow);
    * */

    mp = numroc_( &M, &mb, &myrow, &i_zero, &nprow);
    nq = numroc_( &N, &nb, &mycol, &i_zero, &npcol );
    printf("Thread %d: After mp=%d, np=%d\n",myrank_mpi,mp, nq);

    A_distr = malloc( mp*nq*sizeof(double));
    WORK = (double *)malloc(N*sizeof(double));
    TAU = (double *)malloc(N*sizeof(double));

    // Initialize discriptors (local matrix A is considered as distributed with blocking parameters
    // m, n, i.e. there is only one block - whole matrix A - which is located on process (0,0) )
    lld = MAX( numroc_( &N, &N, &myrow, &i_zero, &nprow ), 1 );
    descinit_( descA, &M, &N, &M, &N, &i_zero, &i_zero, &ictxt, &lld, &INFO );
    lld_distr = MAX( mp, 1 );
    descinit_( descA_distr, &M, &N, &mb, &nb, &i_zero, &i_zero, &ictxt, &lld_distr, &INFO );

    // Call pdgeadd_ to distribute matrix (i.e. copy A into A_distr)
    pdgeadd_( "N", &M, &N, &one, A, &i_one, &i_one, descA, &zero, A_distr, &i_one, &i_one, descA_distr );

    // Now call ScaLAPACK routines
    pdgeqrf_( &M, &N, A_distr, &i_one, &i_one, descA_distr, TAU, WORK, &LWORK, &INFO);

    // Copy result into local matrix
    pdgeadd_( "N", &M, &N, &one, A_distr, &i_one, &i_one, descA_distr, &zero, A, &i_one, &i_one, descA );

    free( A_distr );
    if( myrow==0 && mycol==0 ){
    free( A );
}