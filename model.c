#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <err.h>
#include <lapack.h>
#include "harmonics.h"


extern void dgels(char * trans, int * m, int * n, int * nrhs, double * A, int * lda, double * B, int * ldb, double * work, int * lwork, int * info);

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

int main(int argc, char ** argv)
{
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

	/*************************LAPACK CALL***************************/

	/* define the parameters to give to dgels*/
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