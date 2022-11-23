#ifndef _SCALAPACK_H
#define _SCALAPACK_H


#define BLOCK_CYCLIC_2D 1
#define DLEN_ 9 
#define DTYPE_ 1
#define CTXT_ 2
#define M_ 3 
#define N_ 4 
#define MB_ 5
#define NB_ 6
#define RSRC_ 7 
#define CSRC_ 8
#define LLD_ 9


extern void Cblacs_get(int icontxt, int what, int *val);
extern void Cblacs_pinfo(int * rank, int * nprocs);
extern void Cblacs_setup(int * rank, int * nprocs);
extern void Cblacs_gridinit(int * ictx, const char * order, int nprow, int npcol);
extern void Cblacs_gridinfo(int ictx, int * nprow, int * npcol, int * myrow, int * mycol);
extern void Cblacs_gridexit( int ictx );
extern void Cblacs_exit( int doneflag );

int numroc_(int* N, int* NB, int* IPROC, int *ISRCPROC, int* NPROCS );
int indxg2p_( int* INDXGLOB, int* NB, int* IPROC, int* ISRCPROC, int* NPROCS  );
void pdgemr2d_( int* m, int* n, double* a, int* ia, int* ja, int* desca, double* b, int* ib, int* jb, int* descb, int* ictxt);
void pdgeadd_( char* TRANS, int * M, int * N, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * BETA, double * C, int * IC, int * JC, int * DESCC );
void pdgeqrf_(int* M, int *N, double* A, int* IA, int *JA, int* DESCA, double *TAU, double *WORK, int* LWORK, int *INFO);
void descinit_ (int *desc, const int *m, const int *n, const int *mb, const int *nb, const int *irsrc, const int *icsrc, const int *ictxt, const int *lld, int *info);

void dgeqrf_(int* M, int* N, double *A, int *LDA, double* TAU, double *WORK, int* LWORK, int *INFO);	
void dormqr_(char *side, char *trans, int* M, int* N, int* k, double *A, int *LDA, double* TAU, double* C, int* LDC, double *WORK, int* LWORK, int *INFO);	


#endif