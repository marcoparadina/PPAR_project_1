#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <assert.h>
#include <sys/time.h>
#include <inttypes.h>
#include <string.h>

#include "harmonics.h"

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

/* represent n in <= 8 char  */
void human_format(char * target, long n) {
        if (n < 1000) {
                sprintf(target, "%" PRId64, n);
                return;
        }
        if (n < 1000000) {
                sprintf(target, "%.1f K", n / 1e3);
                return;
        }
        if (n < 1000000000) {
                sprintf(target, "%.1f M", n / 1e6);
                return;
        }
        if (n < 1000000000000ll) {
                sprintf(target, "%.1f G", n / 1e9);
                return;
        }
        if (n < 1000000000000000ll) {
                sprintf(target, "%.1f T", n / 1e12);
                return;
        }
}


void setup_spherical_harmonics(long lmax, struct spherical_harmonics *self)
{
	self->lmax = lmax;
	long sizeCS = (lmax + 1) * (lmax + 1);
	long sizeAB = (lmax + 1) * (lmax + 2) / 2;

	self->CS = malloc(sizeCS * sizeof(double));
	self->A = malloc(sizeAB * sizeof(double));
	self->B = malloc(sizeAB * sizeof(double));
	if (self->CS == NULL || self->A == NULL || self->B == NULL)
		err(1, "cannot allocate harmonics\n");

	/* compute the A, B coefficients */
	for (long l = 2; l <= lmax; l++) {
		double ls = l * l;
		double lm1s = (l - 1) * (l - 1);
		for (long m = 0; m < l - 1; m++) {
			double ms = m * m;
			self->A[PT(l, m)] = sqrt((4 * ls - 1) / (ls - ms));
			self->B[PT(l, m)] = -sqrt((lm1s - ms) / (4 * lm1s - 1));
		}
	}
}

void load_data_points(const char *filename, long npoint, struct data_points *self)
{
	self->npoint = npoint;
	self->phi = malloc(npoint * sizeof(double));
	self->lambda = malloc(npoint * sizeof(double));
	self->V = malloc(npoint * sizeof(double));
	if (self->phi == NULL || self->lambda == NULL || self->V == NULL)
		err(1, "cannot allocate data points\n");

	FILE *f = fopen(filename, "r");
	if (f == NULL)
		err(1, "cannot open %s", filename);

	long nmax=0;
	
	if((strchr(filename,'s')!=NULL) && (strchr(filename,'m')!=NULL) && (strchr(filename,'l')!=NULL) && (strchr(filename,'a')!=NULL)){//small
		nmax=64800;
	}
	else if((strchr(filename,'m')!=NULL) && (strchr(filename,'e')!=NULL) && (strchr(filename,'d')!=NULL) && (strchr(filename,'i')!=NULL) && (strchr(filename,'u')!=NULL)){//medium
		nmax=583200;
	}
	else if((strchr(filename,'h')!=NULL) && (strchr(filename,'i')!=NULL) && (strchr(filename,'g')!=NULL)){//high
		nmax=6480000;
	}
	else if((strchr(filename,'u')!=NULL) && (strchr(filename,'l')!=NULL) && (strchr(filename,'t')!=NULL) && (strchr(filename,'r')!=NULL) && (strchr(filename,'a')!=NULL)){//ultra
			nmax=233280000;
	}

	long step = floor(nmax/npoint);
	step = fmax(step,1);
	//printf("nmax/npoint = step : %d / %d = %d \n",nmax, npoint, step);

	double t1, t2, t3;

	for (long i = 1; i <= step*npoint; i++) {
		if(i%step==0){
			long k = fscanf(f, "%lg %lg %lg", &self->lambda[i/step], &self->phi[i/step], &self->V[i/step]);

			if (k == EOF) {
				if (ferror(f))
					err(1, "read error");
				errx(1, "premature end-of-file after %ld records", i);
			}
			if (k != 3){
				errx(1, "parse error on line %ld", i+1);}
		}
		else{
			long k = fscanf(f,"%lg %lg %lg", &t1, &t2, &t3);

			if (k == EOF) {
				if (ferror(f))
					err(1, "read error");
				errx(1, "premature end-of-file after %ld records", i);
			}
			if (k != 3){

				errx(1, "parse error on line %ld", i+1);}
		}
	}
	fclose(f);
}

void load_spherical_harmonics(const char *filename, long lmax, struct spherical_harmonics *self)
{
	FILE *g = fopen(filename, "r");
	if (g == NULL)
		err(1, "cannot open %s", filename);
	setup_spherical_harmonics(lmax, self);
	for (;;) {
		long l, m;
		double c, s;
		long k = fscanf(g, "%ld %ld %lg %lg", &l, &m, &c, &s);
		if (m == 0 && s != 0)
			errx(1, "non-zero S coefficient with l=%ld and m=0", l);
		self->CS[CT(l, m)] = c;
		if (m > 0)
			self->CS[ST(l, m)] = s;
		if (k == EOF) {
			if (ferror(g))
				err(1, "read error");
			break;
		}
		if (k != 4)
			errx(1, "parse error");
	}
	fclose(g);
}

/*
 * Compute all the (fully normalized) Associated Legendre function of degree <= lmax.
 * On exit, P_{l,m} (with 0 <= m <= l <= lmax) can be found in P[PT(l, m)].
 * P must be preallocated of size (lmax * lmax + 3) / 2.
 * The constants A and B must be precomputed.
 */
void computeP(const struct spherical_harmonics *self, double *P, double sinphi)
{
	double cosphi = sqrt(1 - sinphi * sinphi);
	double temp = 1;
	P[PT(0, 0)] = temp;
	if (self->lmax == 0)
		return;
	P[PT(1, 0)] = sinphi * sqrt(3) * temp;
	temp = 0.5 * sqrt(3) * cosphi * temp;
	P[PT(1, 1)] = temp;
	for (long l = 2; l <= self->lmax; l++) {
		for (long m = 0; m < l - 1; m++)
			P[PT(l, m)] = self->A[PT(l, m)] * (sinphi * P[PT(l - 1, m)] + self->B[PT(l, m)] * P[PT(l - 2, m)]);
		P[PT(l, l - 1)] = sinphi * sqrt(2 * (l - 1) + 3) * temp;
		temp = -sqrt(1.0 + 0.5 / l) * cosphi * temp;
		P[PT(l, l)] = temp;
    //printf("%d",temp);
	}
}

double evaluate(const struct spherical_harmonics *self, const double *P, double lambda)
{
	long lmax = self->lmax;
	long sizeCS = (lmax + 1) * (lmax + 1);
	double scratch[sizeCS];

	for (long l = 0; l <= lmax; l++) {
		/* zonal term */
		scratch[CT(l, 0)] = P[PT(l, 0)];

		/* tesseral terms */
		for (long m = 1; m <= l; m++) {
			scratch[CT(l, m)] = P[PT(l, m)] * cos(m * lambda);
			scratch[ST(l, m)] = P[PT(l, m)] * sin(m * lambda);
		}
	}

	/* dot product */
	double V = 0;
	for (long i = 0; i < sizeCS; i++)
		V += scratch[i] * self->CS[i];
	return V;
}