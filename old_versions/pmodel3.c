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
    MPI_Init(NULL, NULL); //If Cblacs_pinfo is called, MPI_Init is called by it as well if it hasn't been called yet. Since here we never use Cblacs_pinfo,
                          // MPI_Init needs to be called here

    process_command_line_options(argc, argv);

	// Preparations and memory allocation
	int nvar = (lmax + 1) * (lmax + 1);
	//printf("Linear Least Squares with dimension %d x %d\n", npoint, nvar);
	if (nvar > npoint)
		errx(1, "not enough data points");

    //Useful constants
    int i_zero = 0, i_one = 1, i_neg_one = -1;
    double one = 1, zero = 0;
    
    int nprow, npcol, prow, pcol, sysctx;
	
	//Dimensions of process grid
	nprow = NPROW;
	npcol = NPCOL;

	//Cblacs grid initialization
	Cblacs_get(i_zero, i_zero, &sysctx);
	int ctx = sysctx;
	Cblacs_gridinit(&ctx, "Row-Major", nprow, npcol);
	Cblacs_gridinfo(ctx, &nprow, &npcol, &prow, &pcol); // prow, pcol are the row and col of the current process in the process grid

    int M, N, mb_A, nb_A, mb_b, nb_b, mp_A, nq_A, mp_b, nq_b, lld_A, lld_A_distr, lld_b, lld_b_distr;
    int descA[9], descA_distr[9], descb[9], descb_distr[9], lwork, infoA, infoA_distr, infob, infob_distr, info;
    double *A = malloc(npoint*nvar*sizeof(double));
    double *A_distr;
    double *b = malloc(npoint*sizeof(double));
    double *b_distr;
    struct data_points data;

    //Dimensions of global matrix
	M=npoint;
	N=nvar;

	//Blocking parameters (dimention of blocks) respectively for A and b (b is all in one block of dimensions M x 1)
    //POSSIBLE WEAKNESS: do blocking parameters have to be related in some way? Is this a good way to set them?
	mb_A = 10;
	nb_A = 10;
    mb_b = 10;    //POSSIBLE WEAKNESS: We want b to have the same number of rows as A_distr (Is this the way to do it?)
    nb_b = 10;
    
    if((prow == 0) && (pcol == 0)){

        /********************* Build matrices A and b ********************************/	
        //long matrix_size = sizeof(double) * nvar * npoint;
        //char hsize[16];
        //human_format(hsize, matrix_size);
        //printf("Matrix size: %sB\n", hsize);
        /*
        A = malloc(matrix_size);
        if (A == NULL){
            err(1, "cannot allocate matrix");
        } 
        */           

        //printf("Reading data points from %s\n", data_filename);
        
        //Build b
        load_data_points(data_filename, npoint, &data);
        //b = malloc(npoint*sizeof(double));
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

        //Build A
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
        //printf("Built matrices A, b\n");
        //DEBUG
        printf("\nProcess [%d, %d] printing A:\n", prow, pcol);
        /*
        for(int i=0; i<N*M; i++){
            //A[i]=i;
            printf("%d:\t%f\n", i, A[i]);
        }
        */
        for(int i = 0; i<npoint; i++){
            printf("%f\n", b[i]);
        }
        //FLOP = 2. * nvar * nvar * npoint;
        //char hflop[16];
        //human_format(hflop, FLOP);
        //printf("Least Squares (%sFLOP)\n", hflop);
        /***********************************************************************/
	}
    else{
        A = NULL;    
        b = NULL;    
    }
	
	//double start = wtime();  
	
    /**************************** SCALAPACK CALL *********************************/

	//Compute the size of the local part of the matrices A and b, called A_distr and b_distr	
    //POSSIBLE DEFECT: the last arguments might have to be &nprocs, and not &nprows, &npcols. Cf. docs
	mp_A = numroc_(&M, &mb_A, &prow, &i_zero, &nprow);
	nq_A = numroc_(&N, &nb_A, &pcol, &i_zero, &npcol);
    mp_b = mp_A; //numroc_(&M, &mb_b, &prow, &i_zero, &nprow);
    nq_b = numroc_(&i_one, &nb_b, &pcol, &i_zero, &npcol);
    //printf("PG[%d, %d]\tmp_A=%d, nq_A=%d\n ", prow, pcol, mp_A, nq_A);
	A_distr = malloc(mp_A*nq_A*sizeof(double));
	if(A_distr == NULL){
		err(1, "Couldn't allocate A_distr \n");
	}

    b_distr = malloc(mp_b*nq_b*sizeof(double));
    if(b_distr == NULL){
        err(1, "Could not allocate b_distr");
    }	    

	//Initialize the descriptors of A, A_distr, b
	lld_A = M; //MAX( numroc_( &M, &M, &prow, &i_zero, &nprow ), 1 );
    lld_A_distr = MAX(mp_A, 1);
    lld_b = M; //MAX(numroc_(&M, &M, &prow, &i_zero, &nprow), 1); 
    lld_b_distr = MAX(mp_b, 1);
	descinit_(descA, &M, &N, &M, &N, &i_zero, &i_zero, &ctx, &lld_A, &infoA);
	descinit_(descA_distr, &M, &N, &mb_A, &nb_A, &i_zero, &i_zero, &ctx, &lld_A_distr, &infoA_distr);
    descinit_(descb, &M, &i_one, &M, &i_one, &i_zero, &i_zero, &ctx, &lld_b, &infob);
    descinit_(descb_distr, &M, &i_one, &mb_b, &nb_b, &i_zero, &i_zero, &ctx, &lld_b_distr, &infob_distr);
    //printf("Process [%d, %d] initialized its descriptors\n", prow, pcol);
    //DEBUG
    /*
    for(int k = 0; k<9; k++){
        printf("%d: %d\n", k, descb[k]);
    }
    */
    if((infoA != 0) || (infoA_distr != 0) || (infob != 0) || (infob_distr != 0)){
        err(1, "Could not initialize descriptors \n");
    }

    //Distribute matrices A and b
    //POSSIBLE DEFECT: ia, ja, ib, jb might have to be changed, this looks wrong intuitively. 
    int ia, ja, ib, jb;
    ia = i_one;
    ja = i_one;
    ib = i_one;
    jb = i_one;
    //POSSIBLE DEFECT: if nothing else works, probably these are not being copied right. Check the docs about the "context" thing
    //pdgemr2d_(&mp_A, &nq_A, A, &ia, &ja, descA, A_distr, &i_one, &i_one, descA_distr, &ctx); 
    pdgeadd_("N", &M, &N, &one, A, &i_one, &i_one, descA, &zero, A_distr, &i_one, &i_one, descA_distr);
    pdgeadd_("N", &M, &i_one, &one, b, &i_one, &i_one, descb, &zero, b_distr, &i_one, &i_one, descb_distr);
    //printf("Process [%d, %d] has its A_distr\n", prow, pcol);
    //pdgemr2d_(&mp_b, &nq_b, b, &ib, &jb, descb, b_distr, &i_one, &i_one, descb_distr, &ctx);
    //printf("Process [%d, %d] has its b_distr\n", prow, pcol);
    /*
    for(int i = 0; i<mp_b; i++){
        printf("[%d, %d]\t%f\n", prow, pcol, b_distr[i]);
    }
    */
    //DEBUG
    //printf("\nProcess [%d, %d] printing A_distr:\n", prow, pcol);
    
    //DEBUG: Each process prints out its A_distr
    /*
    for(int i=0; i<mp_A*nq_A; i++){
        printf("PG[%d, %d] Before\t%d:\t%f\n", prow, pcol, i, A_distr[i]);
    }
    */    

    //Call the ScaLAPACK least-square-solution routine. The first call is a query to get the optimal value of lwork
    //POSSIBLE DEFECT: ia, ja, ib, jb might have to be changed. Doing this above will change them here as well, this comment is just to notice that those
    //                 parameters are used here as well
    
    double temp_work[1];
    //pdgels_("N", &mp_A, &nq_A, &i_one, A_distr, &ia, &ja, descA_distr, b, &ib, &jb, descb, temp_work, &i_neg_one, &info); //distributing only A
    pdgels_("N", &mp_A, &nq_A, &i_one, A_distr, &ia, &ja, descA_distr, b_distr, &ib, &jb, descb_distr, temp_work, &i_neg_one, &info); //distributing A and b
    //pdgels_("N", &M, &N, &i_one, A, &ia, &ja, descA, b, &ib, &jb, descb, temp_work, &i_neg_one, &info); //no distribution
    if(info!=0){
        err(1, "Query call of pdgels failed \n");
    }
    //printf("Process [%d, %d] finished query PDGELS call\n", prow, pcol);
    lwork = temp_work[0];
    double *work = malloc(lwork * sizeof(double));

    //pdgels_("N", &mp_A, &nq_A, &i_one, A_distr, &ia, &ja, descA_distr, b, &ib, &jb, descb, work, &lwork, &info);  //distributing only A
    pdgels_("N", &mp_A, &nq_A, &i_one, A_distr, &ia, &ja, descA_distr, b_distr, &ib, &jb, descb_distr, work, &lwork, &info);  //distributing A and b
    //pdgels_("N", &M, &N, &i_one, A, &ia, &ja, descA, b, &ib, &jb, descb, work, &lwork, &info); //no distribution
    if(info!=0){
        err(1, "Call of pdgels failed \n");
    }
    //printf("Process [%d, %d] finished actual PDGELS call\n", prow, pcol);
    /*
    for(int i=0; i<mp_A*nq_A; i++){
        printf("PG[%d, %d] After\t%d:\t%f\n", prow, pcol, i, A_distr[i]);
    }
    */
    //Copy the local results into the global matrix
    //POSSIBLE DEFECT: again ia, ja, ib, jb. See above
    //pdgemr2d_(&mp_A, &nq_A, A_distr, &i_one, &i_one, descA_distr, A, &ia, &ja, descA, &ctx); 
    //pdgemr2d_(&mp_b, &nq_b, b_distr, &i_one, &i_one, descb_distr, b, &ib, &jb, descb, &ctx);
    //printf("Process [%d, %d] copied the local result into global matrix\n", prow, pcol);
    pdgeadd_("N", &M, &N, &one, A_distr, &i_one, &i_one, descA_distr, &zero, A, &i_one, &i_one, descA);
    pdgeadd_("N", &M, &i_one, &one, b_distr, &i_one, &i_one, descb_distr, &zero, b, &i_one, &i_one, descb);
    free(A_distr);
    free(b_distr);

    /****************************************************************************/
	
	

    //Save the model in the output file
    if((prow == 0 ) && (pcol == 0 )){
        /*
        //DEBUG:print matrix a
        printf("\nProcess [%d, %d] printing A:\n", prow, pcol);
        for(int i=0; i<N*M; i++){
            //A[i]=i;
            printf("%d:\t%f\n", i, A[i]);
        }
        */        
        /*
        FLOP = 1; //DEBUG
        double t = wtime()  - start;
        double FLOPS = FLOP / t;	
        char hflops[16];
        human_format(hflops, FLOPS);
        //printf("Completed in %.1f s (%s FLOPS)\n", t, hflops);
        */
        for(int i = 0; i< npoint; i++){
            data.V[i] = b[i];
        }
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
        free(b);
    }

    //Exit process grid. This also finalizes MPI	
	Cblacs_gridexit(ctx);
	Cblacs_exit(i_zero);    //This also calls MPI_finalize if the argument is 0
}