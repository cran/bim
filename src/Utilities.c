/*  
 * Utility functions for R/bim package
 *
 * These functions were originally written for QTLCartographer package
 * Modified by Hao Wu <hao@jax.org> August 2003 for the use of R package  
 *
 * The modifications include:
 * 1. Remove all file I/O related functions
 * 2. Use R's function to allocate and free memory
 * 3. Remove all debug related stuff
 *
 */

#include "Main.h"
#include <time.h>

char time_buffer[MAXNAME+1];
char gbuffer[MAXLINE+1];      /*Reusable global buffer space*/
char gname[MAXNAME+1];        /*Reusable global name space*/
char gwarn[MAXNAME+1];
double ranf(long int inix);
long ix, Ix;
int whosemf;
double mapparam;

#define A 16807.0
#define P 2147483647.0
#define a 16807
#define b15 32768
#define b16 65536
#define p 2147483647
#define xnorm 4.656612875E-10

/*
lrtolod converts the likelihood ratio (LR) to the log of the odds (LOD) score.

 LOD = -log10( exp( -LR/2 ) )
 
Because  
  LR  -2ln(L0/L1) 
and 
  LOD = -log10(L0/L1)
*/

double lrtolod(double lr)
{
  double lod;
  lod = - log10( exp( -lr/2.0 ) );  
  return(lod);
}


/* Convert LOD score to LR */
double lodtolr(double lod)
{
  double lr;
  lr  = 2.0*lod* log(10);
  return(lr);
}



/*
  Return a random integer from 1 to in, using ranf(ix) as the
  source of uniform deviates.
  Fortran source from John Monahan,
  and translated into c by Chris Basten 19 January 1994.
*/
long iran(long int xix, long int in)
{
  double rv;
  rv = ranf(xix);
  return ((long) (rv * (double) in) + 1L);
}

/*
 * This subroutine will compute an array of uniformly
 * distributed random numbers.  You need to specify the
 * size of the array, and the numbers will be doubles.
 */

double *ran_arry(int size)
{
  double *arry;
  int ii;
  long I, bound;
  if (size < 0) {
    bound = -size;
    Ix = (long) ((double) P * ranf(bound));
  }
  else
    bound = size;
  arry = dvector(0, bound);
  for (ii = 0; ii <= bound; ii++) {
    I = (Ix * A) / P;
    Ix = (Ix * A - I * P);
    *(arry + ii) = Ix * xnorm;
  }
  return arry;
}

double ranf(long int inix)
{
  long xhi, xalo, leftlo, fhi, k;
  double xx;
  if (inix < 0)
    ix = -1 * inix;
  xhi = ix / b16;
  xalo = (ix - xhi * b16) * a;
  leftlo = xalo / b16;
  fhi = xhi * a + leftlo;
  k = fhi / b15;
  ix = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
  if (ix < 0)
    ix = ix + p;
  xx = (double) ix *xnorm;
  return xx;
}


/*
 * Functions to allocate and free memory
 * Hao Wu removed all error checking stuff, 
 * which will be taken care of by R
 */
void nrerror(char *error_text)
{
  fprintf(stderr, "Numerical Recipes run-time error...\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "...now exiting to system...\n");
  exit(1);
}

float *vector(int nl, int nh)
{
  float *v;
  v = (float *) Calloc((unsigned) (nh - nl + 1), float);
  return v - nl;
}

double *dvector(int nl, int nh)
{
  int ii;
  double *v;

  v = (double *) Calloc( (size_t) (nh - nl + 1) , double);
  return v - nl;
}


int *ivector(int nl, int nh)
{
  int *v,ii;
  v = (int *) Calloc((unsigned) (nh - nl + 1), int);

  return v - nl;
}

char *cvector(int nl, int nh)
{
  char *v;
  int ii;
  
  v = (char *) Calloc((unsigned) (nh - nl + 1), char);

  return v - nl;
}

long *lvector(int nl, int nh)
{
  int ii;
  long *v;
  
  v = (long *) Calloc((unsigned) (nh - nl + 1), long);

  return v - nl;
}

char **cmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i,j;
  char **m;
  
  m = (char **) Calloc((unsigned) (nrh - nrl + 1), char*);
  m -= nrl;
  for (i = nrl; i <= nrh; i++) {
    m[i] = (char *) Calloc((unsigned) (nch - ncl + 1), char);
    m[i] -= ncl;
  }
  /* Do we really need to clear the memory?
    for ( i = nrl; i <= nrh ; i++ )
    for ( j = ncl ; j <= nch ; j++ )
    *(*(m+i)+j) = '\0'; */

  return m;
}


int **imatrix(int nrl, int nrh, int ncl, int nch)
{
  int i, j, **m;
  
  m = (int **) Calloc((unsigned) (nrh - nrl + 1), int*);
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (int *) Calloc((unsigned) (nch - ncl + 1), int);
    m[i] -= ncl;
  }
  /*  for ( i = nrl; i <= nrh ; i++ )
    for ( j = ncl ; j <= nch ; j++ )
    *(*(m+i)+j) = 0; */

  return m;
}

float **matrix(int nrl, int nrh, int ncl, int nch)
{
  int i,j;
  float **m;
  
  m = (float **) Calloc((unsigned) (nrh - nrl + 1), float*);
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (float *) Calloc((unsigned) (nch - ncl + 1), float);
    m[i] -= ncl;
  }
  /* for ( i = nrl; i <= nrh ; i++ )
    for ( j = ncl ; j <= nch ; j++ )
    *(*(m+i)+j) = 0.0; */

  return m;
}



short **smatrix(int nrl, int nrh, int ncl, int nch)
{
  int i,j;
  short **m;
  
  m = (short **) Calloc((unsigned) (nrh - nrl + 1), short*);
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (short *) Calloc((unsigned) (nch - ncl + 1), short);
    m[i] -= ncl;
  }
  /*  for ( i = nrl; i <= nrh ; i++ )
	 for ( j = ncl ; j <= nch ; j++ )
	 *(*(m+i)+j) = 0; */

  return m;
}


void free_vector(float *v, int nl, int nh)
{
  float *tmp;
  tmp = (float *)(v+nl);
  Free(tmp);
}


void free_dvector(double *v, int nl, int nh)
{
  double *tmp;
  tmp = (double *)(v+nl);
  Free(tmp);
}


void free_ivector(int *v, int nl, int nh)
{
  int *tmp;
  tmp = (int *)(v+nl);
  Free(tmp);
}


void free_cvector(char *v, int nl, int nh)
{
  char *tmp;
  tmp = (char *)(v+nl);
  Free(tmp);
}


void free_lvector(long int *v, int nl, int nh)
{
  long int *tmp;
  tmp = (long int *)(v+nl);
  Free(tmp);
}


void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  char *tmp1;
  char **tmp2;
  for (i = nrh; i >= nrl; i--) {
    tmp1 = (char *)(m[i] + ncl);
    Free(tmp1);
  }
  tmp2 = (char **)(m + nrl);
  Free(tmp2);
}


void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  int *tmp1;
  int **tmp2;
  for (i = nrh; i >= nrl; i--) {
    tmp1 = (int *)(m[i] + ncl);
    Free(tmp1);
  }
  tmp2 = (int **)(m + nrl);
  Free(tmp2);
}


void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  float *tmp1;
  float **tmp2;
  for (i = nrh; i >= nrl; i--) {
    tmp1 = (float *)(m[i] + ncl);
    Free(tmp1);
  }
  tmp2 = (float **)(m + nrl);
  Free(tmp2);
}


void free_smatrix(short int **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  short *tmp1;
  short **tmp2;
  for (i = nrh; i >= nrl; i--) {
    tmp1 = (short *)(m[i] + ncl);
    Free(tmp1);
  }
  tmp2 = (short **)(m + nrl);
  Free(tmp2);
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i,j;
  double **m;
  
  m = (double **) Calloc((unsigned) (nrh - nrl + 1), double*);
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = (double *) Calloc((unsigned) (nch - ncl + 1), double);
    m[i] -= ncl;
  }

  return m;
}


void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  double *tmp1;
  double **tmp2;
  for (i = nrh; i >= nrl; i--) {
    tmp1 = (double *)(m[i] + ncl);
    Free(tmp1);
  }
  tmp2 = (double **)(m + nrl);
  Free(tmp2);
}


/*****************************************************************************
 *        This is a general Gamma variate generator with the shape
 *     parameter beta > 0.25. Coded from the Algorithm GBH of Cheng and Feast
 *     (1980 Communications of the ACM 23:389-394) by Zhao-Bang Zeng on Jan. 30,
 *     1992.
 *     This was translated into c by Chris Basten, 19 Jan. 1994.
 *****************************************************************************/
double gamgbh(double beta, long int xix)
{
  double aa, bb, c, d, t, h1, h2, u, u1, u2, w, tmp;

  if (beta <= 0.25)
    nrerror("Error in gamgbh:  beta <= 0.25...");
  aa = beta - 0.25;
  bb = beta / aa;
  c = 2.0 / aa;
  d = c + 2.0;
  t = 1.0 / sqrt(beta);
  h1 = (0.4417 + 0.0245 * t / beta) * t;
  h2 = (0.222 - 0.043 * t) * t;
  do {
    u = ranf(xix);
    u1 = ranf(xix);
    u2 = u1 + h1 * u - h2;
    if (u2 <= 0.0 || u2 >= 1.0)
      tmp = 1.0;
    else {
      w = bb * pow((u1 / u2), 4.0);
      if (c * u2 - d + w + 1.0 / w <= 0.0)
	return (aa * w);
      else
	tmp = c * log(u2) - log(w) + w - 1.0;
    }
  } while (tmp >= 0.0);
  return (aa * w);
}


/****************************************************************************
      This gamma generator generates gamma variates by composition-
   rejection from the Weibull density for the shape parameter beta
   smaller than 1 (0 < beta < 1). Coded from I. Vaduva (1977 Math.
   Operationsforsch. Statist., Ser. Statistics, Vol. 8:545-576) by
   Zhao-Bang Zeng on Feb. 5 1992.   (Best one of this kind)

   Translated into c by Chris Basten on 19 January 1994.
 *****************************************************************************/
double gamnl1(double beta, double aa, double bb, double pp, long int xix)
{
  double gamnl1r, u1, u2, u0;

  if (beta >= 1.0 || beta <= 0.0)
    nrerror("beta in gamnl1 is not in (0,1)");

/*  Change > to < in the following line suggested by Ian Painter. 21 Oct. 1996*/
  if (ranf(xix) < pp) {
    u0 = ranf(xix);
    gamnl1r = u0 = pow(u0, aa);
    u1 = ranf(xix);
    while (u0 >= u1) {
      u2 = ranf(xix);
      if (u1 < u2) {
	u0 = ranf(xix);
	gamnl1r = u0 = pow(u0, aa);
      }
      else
	u0 = u2;
      u1 = ranf(xix);
    }
  }
  else
    do {
      u0 = ranf(xix);
      u0 = pow(u0, bb);
      gamnl1r = 1.0 - log(ranf(xix));
    } while (gamnl1r >= u0);
  return (gamnl1r);
}


/*
For:
  MIN_INT < lb < ub < MAX_INT  and ub - lb < MAX_INT
Do:
  Shuffle  v = ivector(lb,ub)
Where:
  ranf() gives a uniform random number in [0.0, 1.0)

foreach i in [lb,ub) in order, pick a random j from (i,ub]
and switch v(i), v(j)

*/
void shuffle_ivector(int *v, int lb, int ub)
{
  int i, j, k;
  for (i = lb; i < ub; i++) {
    j = i + (int) ((double) (ub - i) * ranf(i)) + 1;
    k = *(v + i);
    *(v + i) = *(v + j);
    *(v + j) = k;
    if (j > ub || j <= i)
      printf("\nlb = %d, i = %d, j = %d, ub = %d\n", lb, i, j, ub);
  }
}


int get_int(void)
{
  char buffer[15], ch;
  int ii, ans;
  for (ii = 0; ii < 15; ii++)
    *(buffer + ii) = '\0';
  while (isspace(ch = getchar()));
  *(buffer + 0) = ch;
  for (ii = 1; ii < 15 && *(buffer + ii - 1) != '\n'; ii++)
    *(buffer + ii) = getchar();
  if (ii < 15) {
    *(buffer + ii - 1) = '\0';
    ans = atoi(buffer);
  }
  else
    ans = -1;
  return (ans);
}

/*
This is a function to convert between recombination
frequencies and distance in Morgans or centiMorgans.

value is what will be converted.
flag indicates how it will be converted:
flag =
   -2  => value cM to rvalue Rec. Freq.
   -1  => value  M to rvalue Rec. Freq.
    0  => value is returned unchanged.
    1  => value Rec. Freq. to rvalue M
    2  => value Rec. Freq. to rvalue cM


The global variable whosemf determines which mapfunction to
use.

  whosemf =
             1 => Haldane (1919)  
             2 => Kosambi (1944)  
             3 => Morgan (1928) (Fixed)  
             4 => Carter and Falconer (1951)
             5 => Rao et al (1977)
             6 => Sturt (1976)
             7 => Felsenstein (1979)
             8 => Karlin (1984)
See Ben Hui Liu (1998) "Statistical Genomics: Linkage, Mapping and QTL Analysis" p319
CRC Press

Send Morgans or r to the following functions:


Kosambi, iKosambi
Haldane, iHaldane
Morgan, iMorgan
CarterFalconer, iCarterFalconer
Rao, iRao
Sturt, iSturt
Felsenstein, iFelsenstein
Karlin, iKarlin


Check that 
  0.0 < m 
  0.0 < r < 0.5
*/
double mapfunc(double value, int flag)
{
  double rvalue;
  double mval,rval;
  /* 0.0 < value
     0.0 < value < 0.5 if value is a recombination frequency */
  if ( value < 0.0 )
    return(-1.0);  /* -1 for a negative distance/probability */
  else if ( value == 0.0 ) 
    return(0.0);    
  if ( flag > 0 && value >= 0.5 ) /* -2 for a rec. prob. >= 1/2 */
    return(-2.0);
  
    
  mval = rval = value;
  if (flag == -2 )
    mval = mval*0.01;
  else if ( flag == 0)
    return (value);
   
  if ( flag < 0 ) {
    switch(whosemf) {
      default: case 1: rval = iHaldane(mval);  break;
      case 2: rval = iKosambi(mval);  break;
      case 3: rval = iMorgan(mval);  break;
      case 4: rval = iCarterFalconer(mval);  break;
      case 5: rval = iRao(mval);  break;
      case 6: rval = iSturt(mval);  break;
      case 7: rval = iFelsenstein(mval);  break;
      case 8: rval = iKarlin(mval);  break;
    }
  
  }
  else if ( flag > 0 ) {
    switch(whosemf) {
      default: case 1: mval = Haldane(rval);  break;
      case 2: mval = Kosambi(rval);  break;
      case 3: mval = Morgan(rval);  break;
      case 4: mval = CarterFalconer(rval);  break;
      case 5: mval = Rao(rval);  break;
      case 6: mval = Sturt(rval);  break;
      case 7: mval = Felsenstein(rval);  break;
      case 8: mval = Karlin(rval);  break;
    }  
  }
  
    
  if (flag == 2 )  /*change Morgans to Centimorgans*/
    rvalue = mval*100.0;
  else if ( flag == 1 )
    rvalue = mval;
  else if ( flag == -1 || flag == -2)
    rvalue = rval;
  else 
    rvalue = value;
  
  return (rvalue);
}

/*Inverse of Kosambi mapping function*/
double iKosambi(double mm)
{
  double rr;   
  rr = 0.5 * (1.0 - exp(-4.0 * mm)) / (1.0 + exp(-4.0 * mm));
  return (rr);
}
/*Kosambi mapping function*/
double Kosambi(double rr)
{
  double mm;   
  mm = 0.25 * log((1.0 + 2.0 * rr) / (1.0 - 2.0 * rr));
  return (mm);
}

/*Inverse of Morgan mapping function*/
double iMorgan(double mm)
{
  double rr;   
  rr = mm;
  return (rr);
}
/*Morgan mapping function*/
double Morgan(double rr)
{
  double mm;   
  mm = rr;
  return (mm);
}

/*Inverse of Hadane mapping function*/
double iHaldane(double mm)
{
  double rr;   
  rr = 0.5 * (1.0 - exp(-2.0 * mm));;
  return (rr);
}
/*Hadane mapping function*/
double Haldane(double rr)
{
  double mm;   
  mm = -0.5 * log(1.0 - 2.0 * rr);
  return (mm);
}
/*Inverse of CarterFalconer mapping function*/
double iCarterFalconer(double mm)
{
  double rr,ru,rl,delta,mt;   
      delta = (double) MAPDELTA;
      rl = delta;
      ru = 0.5-delta;
      do {
        rr = 0.5*(rl+ru);
        mt = CarterFalconer(rr);
        if ( mt > mm )
          ru = rr;
        else
          rl = rr;
      } while ( fabs(mt-mm) > delta );
  return (rr);
}
/*CarterFalconer mapping function*/
double CarterFalconer(double rr)
{
  double mm;    
  mm = 0.5*( atan(2.0*rr) + 0.5*log( (1.0+2.0*rr)/(1.0-2.0*rr) ) ) ;
  return (mm);
}
/*Inverse of Felsenstein mapping function*/
double iFelsenstein(double mm)
{
  double rr,kk;   
  kk = mapparam;   
  rr = (1.0-exp(2.0*(kk-2.0)*mm))/(2.0*(1.0-(kk-1.0)*exp(2.0*(kk-2.0)*mm)));
  return (rr);
}
/*Felsenstein mapping function*/
double Felsenstein(double rr)
{
  double mm,kk;
  kk = mapparam;   /* kk can't be 2 */
  mm = log( (1.0-2.0*rr)/(1.0-2.0*(kk-1.0)*rr) )/(2.0*(kk-2.0)) ;
  return (mm);
}
/*Inverse of Karlin mapping function*/
double iKarlin(double mm)
{
  double rr,nn;
  nn = mapparam;   
  rr = 0.5*(1- pow((1.0-2.0*mm/nn),nn) );
  return (rr);
}
/*Karlin mapping function*/
double Karlin(double rr)
{
  double mm,nn;
  nn = mapparam;      
  mm = 0.5*nn*(1.0-pow((1.0-2.0*rr),1.0/nn));
  return (mm);
}
/*Inverse of Rao mapping function*/
double iRao(double mm)
{
  double rr,ru,rl,delta,mt;   
      delta = (double) MAPDELTA;
      rl = delta;
      ru = 0.5-delta;
      do {
        rr = 0.5*(rl+ru);
        mt = Rao(rr);
        if ( mt > mm )
          ru = rr;
        else
          rl = rr;
      } while ( fabs(mt-mm) > delta );
  return (rr);
}
/*Rao mapping function*/
double Rao(double rr)
{
  double mm,pp;
  pp = mapparam;     
  mm = pp*(2.0*pp-1.0)*(1.0-4.0*pp)*log(1.0-2.0*rr)/6.0 + (8.0*pp*(pp-1.0)*(2.0*pp-1.0)*atan(2.0*rr) +  pp*(1.0-pp)*(4.0*pp+1.0) *log((1.0+2.0*rr)/(1.0-2.0*rr)))/3.0 + (1.0-pp)*(1.0-2.0*pp)*(1.0-4.0*pp)*rr;
  return (mm);
}

/*Inverse of Sturt mapping function*/
double iSturt(double mm)
{
  double rr,ll;
  ll = (double) mapparam;   
  rr = 0.5*(1.0-(1.0-mm/ll)*exp(mm*(1.0-2.0*ll)/ll) );
  return (rr);
}

/*Sturt mapping function*/
double Sturt(double rr)
{
  double mm,ml,mu,rt,delta;   
  ml = delta = (double) MAPDELTA;
  do {
    mu = ml+1.0;
    rt = iSturt(mu);
    if (rt < rr )
      ml = ml + 1.0;
  
  } while ( rt < rr );
  do {
        mm = 0.5*(ml+mu);
        rt = iSturt(mm);
        if ( rt > rr )
          mu = mm;
        else
          ml = mm;
  } while ( fabs(rt-rr) > delta );
  return (mm);
}


long get_a_seed(void)
{
  time_t tptr;
  time(&tptr);
#if defined(THINK_C)
  tptr = tptr / 2L;
#endif
  return ((long) tptr);
}


char *asctime2(void)
{
/*  static char time_buffer[MAXNAME];*/
  int i;
  time_t tptr;
  struct tm *tms;
  size_t len;
  for ( i=0;i<MAXNAME;i++ )
    time_buffer[i] = '\0';
  time(&tptr);
  tms = localtime(&tptr);
  len = strftime(time_buffer, MAXNAME, "%H:%M:%S on %A, %d %B %Y\n", tms);
#if defined(THINK_C) 
    return time_buffer;
#else
  if (len == 0)
    return NULL;
  else
    return time_buffer;
#endif
}


#if defined(DSIGN)
double dsign(double val1, double val2)
{
  double xx;
  xx = fabs(val1);
  if (val2 < 0.0)
    xx = -xx;
  return (xx);
}
#endif



/*
 * dtranspose does an arbitrary transpose of one matrix
 * onto another. You must have allocated the memory that
 * mm1 and mm2 have pointed to.  lr, and lc are the initial
 * row and column while ur and uc are the final row and
 * column for the transpositon.  Note that mm1 =
 * dmatrix(lr,ur,lc,uc); mm2 = dmatrix(lc,uc,lr,ur); at
 * least.  These can be bigger, so that one can transpose
 * an arbitrary part of one matrix onto another.
 *
 * By Chris Basten, January 1994
 */
int dtranspose(double **mm1, double **mm2, int lr, int lc, int ur, int uc)
{
  int ii, jj;
  if (lr >= ur || lc >= uc)
    return (1);
  for (ii = lr; ii <= ur; ii++)
    for (jj = lc; jj <= uc; jj++)
      *(*(mm2 + jj) + ii) = *(*(mm1 + ii) + jj);
  return (0);
}



#if defined(ITOA)


/*
 * I have found some of the functions on other machines,
 * namely the Decstation.  Also, the GNU C compiler has
 * them.
 */
char *
     strupr(char *s)
{
  int len, ii;
  if ( s == NULL )
    return(s);
  len = strlen(s);
  for (ii = 0; ii < len; ii++)
    if ( islower(s[ii]) )
      s[ii] = toupper(s[ii]);
  return s;
}

char *strlwr(char *s)
{
  int len, ii;
  if ( s == NULL )
    return(s);
  len = strlen(s);
  for (ii = 0; ii < len; ii++) 
    if ( isupper(s[ii]) )
      s[ii] = tolower(s[ii]);
  return(s);
}

#endif


#if defined(DIVT)

div_t div(int numer, int denom);
{
  div_t x;
  x.quot = 0;
  x.rem = 0;
  if (denom != 0) {
    x.quot = numer / denom;
    x.rem = numer - x.quot * denom;
  }
  else
    fprintf(stderr, "\n\nSorry, divide by zero in _div...\n");
  return x;
}

ldiv_t ldiv(long numer,long denom);
{
  ldiv_t x;
  x.quot = 0;
  x.rem = 0;
  if (denom != 0) {
    x.quot = numer / denom;
    x.rem = numer - x.quot * denom;
  }
  else
    fprintf(stderr, "\n\nSorry, divide by zero in _ldiv...\n");
  return x;
}
#endif

