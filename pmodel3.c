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

int lmax = -1;
int npoint;
char * data_filename;
char * model_filename;


int MAX(int a, int b){
	if(a > b){
		return a;
	}
	return b;
}

int MOD(int a, int b){
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
	/* missing required args? */
	if (data_filename == NULL || model_filename == NULL || lmax < 0 || npoint <= 0)
		usage(argv);
}


//TODO: b needs to be distributed as well, otherwise least-square can't be computed because it would have incompatible dimensions between
//      the A_distr matrix and b. The least-square-solution, found by pdgels, minimizes |A_distr*x=b|. To make dimensions compatible b needs
//      to be replaced with a b_distr
int main(int argc, char ** argv){

    MPI_Init(NULL, NULL);

    process_command_line_options(argc, argv);

	// Preparations and memory allocation
	int nvar = (lmax + 1) * (lmax + 1);
	//printf("Linear Least Squares with dimension %d x %d\n", npoint, nvar);
	if (nvar > npoint)
		errx(1, "not enough data points");

    //Useful constants
    int i_zero = 0, i_one = 1, i_neg_one = -1;
    
    int nprow, npcol, prow, pcol, sysctx;
	

	//Dimensions of process grid
	nprow = NPROW;
	npcol = NPCOL;

	//Cblacs grid initialization
	Cblacs_get(i_zero, i_zero, &sysctx);
	int ctx = sysctx;
	Cblacs_gridinit(&ctx, "Row-Major", nprow, npcol);
	Cblacs_gridinfo(ctx, &nprow, &npcol, &prow, &pcol); // prow, pcol are the row and col of the current process in the process grid

    int M, N, mb_A, nb_A, mb_b, nb_b, mp, nq, lld, lld_distr, info, descA[9], descA_distr[9], descb[9], lwork, infoA, infoA_distr, infob, lld_b;
    double FLOP;
    double *A;
    double *A_distr;
    double *b;
    struct data_points data;

    //Dimensions of global matrix
	M=npoint;
	N=nvar;

	//Blocking factors (dimention of blocks) resp. for A and b (b is all in one block of dimensions M x 1)
	mb_A = 10;
	nb_A = 10;
    mb_b = M;
    nb_b = 1;
    
    if((prow == 0) && (pcol == 0)){

        /********************* Build matrix A ********************************/	
        long matrix_size = sizeof(double) * nvar * npoint;
        char hsize[16];
        human_format(hsize, matrix_size);
        //printf("Matrix size: %sB\n", hsize);

        A = malloc(matrix_size);
        if (A == NULL)
            err(1, "cannot allocate matrix");

        //printf("Reading data points from %s\n", data_filename);
        
        load_data_points(data_filename, npoint, &data);
        b = malloc(npoint*sizeof(double));
        for(int i = 0; i< npoint; i++){
            b[i] = data.V[i];
        }
        //printf("Successfully read %d data points\n", npoint);


        double * P = malloc((lmax + 1) * (lmax + 1) * sizeof(*P));
        double * v = malloc(npoint * sizeof(*v));
        if (P == NULL || v == NULL)
            err(1, "cannot allocate data points\n");

        //printf("Building matrix\n");
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
        FLOP = 2. * nvar * nvar * npoint;
        char hflop[16];
        human_format(hflop, FLOP);
        //printf("Least Squares (%sFLOP)\n", hflop);
        /***********************************************************************/
	}
    else{
        A = NULL;    
        b = NULL;    
    }
	
	double start = wtime();  
	
    /**************************** SCALAPACK CALL *********************************/

	//Compute the size of the local part of the matrix A, called A_distr	
	mp = numroc_(&M, &mb_A, &prow, &i_zero, &nprow);
	nq = numroc_(&N, &nb_A, &pcol, &i_zero, &npcol);

	A_distr = malloc(mp*nq*sizeof(double));
	if(A_distr == NULL){
		err(1, "Couldn't allocate A_distr \n");
	}	

	//Initialize the descriptors of A, A_distr, b
	lld = MAX( numroc_( &M, &M, &prow, &i_zero, &nprow ), 1 );
    lld_distr = MAX(numroc_(&M, &mb_A, &prow, &i_zero, &nprow), 1);
    lld_b = MAX(numroc_(&M, &mb_b, &prow, &i_zero, &nprow), 1); 
	descinit_(descA, &M, &N, &M, &N, &i_zero, &i_zero, &ctx, &lld, &infoA);
	descinit_(descA_distr, &M, &N, &mb_A, &nb_A, &i_zero, &i_zero, &ctx, &lld_distr, &infoA_distr);
    descinit_(descb, &M, &i_one, &mb_b, &nb_b, &i_zero, &i_zero, &ctx, &lld_b, &infob);
    //DEBUG
    /*
    for(int k = 0; k<9; k++){
        printf("%d: %d\n", k, descb[k]);
    }
    */
    if((infoA != 0) || (infoA_distr != 0) || (infob != 0)){
        err(1, "Could not initialize descriptors \n");
    }

    //Distribute the matrix
    //TODO: ia, ja might have to be changed, this looks wrong intuitively
    int ia, ja, ib, jb;
    ia = i_one;
    ja = i_one;
    ib = i_one;
    jb = i_one;
    pdgemr2d_(&mp, &nq, A, &ia, &ja, descA, A_distr, &i_one, &i_one, descA_distr, &ctx); //parameters 3,4 don't really make sense to me but they seem to work

    //Call the ScaLAPACK least-square-solution routine. The first call is a query to get the optimal value of lwork
    double temp_work[1];
    pdgels_("N", &mp, &nq, &i_one, A_distr, &ia, &ja, descA_distr, b, &ib, &jb, descb, temp_work, &i_neg_one, &info);
    printf("INFO: %d", info);
    if(info!=0){
        err(1, "Query call of pdgels failed \n");
    }
    lwork = temp_work[0];
    double *work = malloc(lwork * sizeof(double));

    pdgels_("N", &mp, &nq, &i_one, A_distr, &ia, &ja, descA_distr, b, &ib, &jb, descb, work, &i_neg_one, &info);
    if(info!=0){
        err(1, "Call of pdgels failed \n");
    }

    //Copy the local results into the global matrix
    pdgemr2d_(&mp, &nq, A_distr, &i_one, &i_one, descA_distr, A, &ia, &ja, descA, &ctx); //parameters 7,8 don't really make sense to me but they seem to work

    free(A_distr);

    /****************************************************************************/
	
	

    //Save the model in the output file
    if((prow == 0 ) && (pcol == 0 )){
        FLOP = 1; //DEBUG
        double t = wtime()  - start;
        double FLOPS = FLOP / t;	
        char hflops[16];
        human_format(hflops, FLOPS);
        //printf("Completed in %.1f s (%s FLOPS)\n", t, hflops);
        double res = 0;
        for (int j = nvar; j < npoint; j++){
            res += data.V[j] * data.V[j];
        }		
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
        free(A);
    }

    //Exit process grid. This also finalizes MPI	
	Cblacs_gridexit(ctx);
	Cblacs_exit(i_zero);
    //MPI_Finalize();	
}