#ifndef MYLINPACK_GUARD
#define MYLINPACK_GUARD
/*
 * The subroutines below were originally written in
 * Fortran. They appeared in Dongarra, J.J, Moler, C.B.,
 * Bunch, J.R. and Stewart G.W. 1979 LINPACK Users' Guide,
 * SIAM, Philadelphia.
 *
 * The C translations were done by Christopher J. Basten at
 * North Carolina State University in December of 1993.
 * These are simple ports in that they do not unroll loops
 * when the increments are both one.
 *
 */




int sqrst(double **xx, int ldx, int nn, int pp, double *yy, double tol, double *bb, double *rsd, int *kk, int *jpvt, double *qraux);
int sqrdc(double **yy, int ldx, int nn, int pp, double *qraux, int *jpvt, int job);
int sqrsl(double **yy, int ldx, int nn, int pp, int kk, double *qraux, double *y, double *qy, double *qty, double *bb, double *rsd, double *xb, long int job);
int spodi(double **yy, int lda, int n, double *det, int job);
int strsl(double **yy, int ldt, int n, double *b, int job);

void dcswap(double *a, double *b);
int ssidi(double **a, int lda, int n, int *kpvt, double *det, int *inert, double *work, int job);
int ssifa(double **a, int lda, int n, int *kpvt);
int ssico(double **a, int lda, int n, int *kpvt, double *rcond, double *z);
int spofa(double **a, int lda, int n);

double amax1(double x, double y);

#endif

