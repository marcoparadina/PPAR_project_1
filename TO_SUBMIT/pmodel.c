#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include "harmonics.h"
#include <mpi.h>
#include "scalapack.h"

//Dimensions of process grid
#define NPROW 2
#define NPCOL 2

long lmax = -1;
long npoint;
char * data_filename;
char * model_filename;


long max(long a, long b){
	if(a > b){
		return a;
	}
	return b;
}

long min(long a, long b){
	if(a < b){
		return a;
	}
	return b;
}

long MOD(long a, long b){
	return (a%b);	
}

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
	
	if (data_filename == NULL || model_filename == NULL || lmax < 0 || npoint <= 0)
		usage(argv);
}

int main(int argc, char ** argv){
    
    //If Cblacs_pinfo is called, MPI_Init is called by it as well if it hasn't been called yet. Since here we never use Cblacs_pinfo,
    // MPI_Init needs to be called here
    MPI_Init(NULL, NULL); 

    process_command_line_options(argc, argv);

	//Preparations
	long nvar = (lmax + 1) * (lmax + 1);
	if (nvar > npoint)
		errx(1, "not enough data points");

    //Useful constants
    long i_zero = 0, i_one = 1;
    long i_neg_one = -1;
    
    int nprow, npcol, prow, pcol, sysctx;
	
	//Dimensions of process grid
	nprow = NPROW;
	npcol = NPCOL;

	//Cblacs grid initialization: prow, pcol are the row and col of the current process in the process grid
	Cblacs_get(i_zero, i_zero, &sysctx);
	int ctx = sysctx;
	Cblacs_gridinit(&ctx, "Row-Major", nprow, npcol);
	Cblacs_gridinfo(ctx, &nprow, &npcol, &prow, &pcol);

    //Declarations and memory allocations
    long M, N, mb_A, nb_A, mb_b, nb_b, mp_A, nq_A, mp_b, nq_b, lld_A, lld_A_distr, lld_b, lld_b_distr, lwork;
    int descA[9], descA_distr[9], descb[9], descb_distr[9], infoA, infoA_distr, infob, infob_distr, info;
    double start = 0, t;
    double *A;
    double *A_distr;
    double *b;
    double *b_distr;
    struct data_points data;

    //Dimensions of global matrix
	M=npoint;
	N=nvar;

	//Blocking parameters (dimention of blocks) respectively for A and b
	mb_A = 50;
	nb_A = 50;
    mb_b = 50; 
    nb_b = 50;
    
    if((prow == 0) && (pcol == 0)){

        /********************* Build matrices A and b ********************************/	
        
        //Build b
        load_data_points(data_filename, npoint, &data);
        b = malloc(npoint*sizeof(double));

        //Copy data.V into b
        for(long i = 0; i< npoint; i++){
            b[i] = data.V[i];
        }
        
        double * P = malloc((lmax + 1) * (lmax + 1) * sizeof(*P));
        double * v = malloc(npoint * sizeof(*v));
        if (P == NULL || v == NULL)
            err(1, "cannot allocate data points\n");

        //printf("Building matrix\n");
        struct spherical_harmonics model;
        setup_spherical_harmonics(lmax, &model);

        //Build A
        A = malloc(M*N*sizeof(double));
        for (long i = 0; i < npoint; i++) {
            computeP(&model, P, sin(data.phi[i]));
            
            for (long l = 0; l <= lmax; l++) {
                /* zonal term */
                A[i + npoint * CT(l, 0)] = P[PT(l, 0)];
        
                /* tesseral terms */
                for (long m = 1; m <= l; m++) {
                    A[i + npoint * CT(l, m)] = P[PT(l, m)] * cos(m * data.lambda[i]);
                    A[i + npoint * ST(l, m)] = P[PT(l, m)] * sin(m * data.lambda[i]);
                }
            }
        }

        /***********************************************************************/
        start = wtime();
	}
    else{
        A = NULL;    
        b = NULL;    
    }
	
    /**************************** SCALAPACK CALLS *********************************/

	//Compute the size of the local part of the matrices A and b, called A_distr and b_distr	
	mp_A = numroc_(&M, &mb_A, &prow, &i_zero, &nprow);
	nq_A = numroc_(&N, &nb_A, &pcol, &i_zero, &npcol);
    mp_b = numroc_(&M, &mb_b, &prow, &i_zero, &nprow);
    nq_b = numroc_(&i_one, &nb_b, &pcol, &i_zero, &npcol);  
    
	A_distr = malloc(mp_A*nq_A*sizeof(double));
	if(A_distr == NULL){
		err(1, "Couldn't allocate A_distr \n");
	}

    b_distr = malloc(mp_b*nq_b*sizeof(double));
    if(b_distr == NULL){
        err(1, "Could not allocate b_distr");
    }	    

	//Initialize the descriptors of A, A_distr, b, b_distr
	lld_A = M;
    lld_A_distr = max(mp_A, 1);
    lld_b = M;
    lld_b_distr = max(mp_b, 1);
	descinit_(descA, &M, &N, &M, &N, &i_zero, &i_zero, &ctx, &lld_A, &infoA);
	descinit_(descA_distr, &M, &N, &mb_A, &nb_A, &i_zero, &i_zero, &ctx, &lld_A_distr, &infoA_distr);
    descinit_(descb, &M, &i_one, &M, &i_one, &i_zero, &i_zero, &ctx, &lld_b, &infob);
    descinit_(descb_distr, &M, &i_one, &mb_b, &nb_b, &i_zero, &i_zero, &ctx, &lld_b_distr, &infob_distr);
    if((infoA != 0) || (infoA_distr != 0) || (infob != 0) || (infob_distr != 0)){
        err(1, "Could not initialize descriptors \n");
    }

    //Distribute matrices A and b into local matrices A_distr, b_distr
    int ia, ja, ib, jb;
    ia = i_one;
    ja = i_one; 
    ib = i_one;
    jb = i_one;
   
    Cpdgemr2d(M,N, A, ia, ja, descA, A_distr, i_one, i_one, descA_distr, ctx);  
    Cpdgemr2d(M,i_one, b, ib, jb, descb, b_distr, i_one, i_one, descb_distr, ctx);        
    
    //Call to the ScaLAPACK least-square-solution routine. The first call is a query to get the optimal value of lwork.
    double temp_work[1];
    pdgels_("N", &M, &N, &i_one, A_distr, &i_one, &i_one, descA_distr, b_distr, &i_one, &i_one, descb_distr, temp_work, &i_neg_one, &info);
    if(info!=0){
        err(1, "Query call of pdgels failed \n");
    }
    lwork = temp_work[0];
    double *work = malloc(lwork * sizeof(double));

    pdgels_("N", &M, &N, &i_one, A_distr, &i_one, &i_one, descA_distr, b_distr, &i_one, &i_one, descb_distr, work, &lwork, &info);
    if(info!=0){
        err(1, "Call of pdgels failed \n");
    }
    

    //Copy the least-square solution into the global matrix b (A doesn't need to be updated: the LSS is contained in b)
    Cpdgemr2d(M,i_one, b_distr, i_one, i_one, descb_distr, b, ib, jb, descb, ctx); 
    free(A_distr);
    free(b_distr);

    /****************************************************************************/	

    //Save the model in the output file
    if((prow == 0 ) && (pcol == 0 )){

        //The first process prints the time that it took to solve the least-square-problem. The first process is the one that takes the longest, 
        //since it has to handle the global matrices, so its execution time represents the execution time of the whole program
        t = wtime() - start;
        printf("Completed in %.1f s\n", t);

        //Copy b back into data.V
        for(long i = 0; i< npoint; i++){
            data.V[i] = b[i];
        }
        double res = 0;
        for (long j = nvar; j < npoint; j++){
            res += data.V[j] * data.V[j];
        }		
        printf("residual sum of squares %g\n", res);

        printf("Saving model in %s\n", model_filename);
        FILE *g = fopen(model_filename, "w");
        if (g == NULL)
            err(1, "cannot open %s for writing\n", model_filename);
        for (long l = 0; l <= lmax; l++) {
            fprintf(g, "%ld\t0\t%.18g\t0\n", l, data.V[CT(l, 0)]);
            for (long m = 1; m <= l; m++)
                fprintf(g, "%ld\t%ld\t%.18g\t%.18g\n", l, m, data.V[CT(l, m)], data.V[ST(l, m)]);
        }

        //Prints the execution time and residual sum of squares into timer_file
        FILE *e = fopen("timer_file", "w");
        if (e == NULL)
            err(1, "cannot open %s for writing\n", "timer_file");
        fprintf(e, "pmodel3_clean completed its execution in %fs\nResidual sum of squares %g", t, res);
            
        free(A);
        free(b);
    }

    //Exit process grid. This also finalizes MPI	
	Cblacs_gridexit(ctx);
	Cblacs_exit(i_zero); 
}