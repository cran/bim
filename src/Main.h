#ifndef MAIN_GUARD
#define MAIN_GUARD
/*  
          Main.h 
          
Copyright (C) 1996 Christopher J. Basten, Bruce S. Weir and Zhao-Bang Zeng.

This file is part of QTL Cartographer. QTL Cartographer is free software; you
can redistribute it and/or modify it under the terms of the GNU  General
Public License as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

QTL Cartographer is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
QTL Cartographer; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include "LocalD.h"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include <R.h>

#if defined (__WIN32__)
#include <conio.h>
#endif

#if defined(BORLAND)
#include <windows.h>
void WYield();
#endif

#if defined(THINK_C)   
#include <sioux.h>
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#endif


/*  These are some default parameters for Simulations...*/
#define MAX_REPS 10001      /*By default, only up to 10,000 repetitions for permutations or bootstraps. Probably shouldn't go higher than 31,000. */
#define REPS 0              /*No reps as starting point Various*/
#define NN 200              /*Sample size   Rcross */
#define NINBRED  1          /*Vestige of something ? */
#define DEFsigl  0.0        /*StDev of markers/chromosome Rmap*/
#define DEFsigs  0.0        /*StDev of intermarker distances Rmap*/
#define DEFm     4          /*Number of chromosomes Rmap*/
#define DEFqtl   9          /*Ave # of QTL per trait Rqtl*/
#define DEFl     16         /*Ave # of markers/chromosome Rmap*/
#define DEFs     10.0       /*Ave intermarker distance   Rmap*/
#define DEFbrdrs 0.0        /*Ave amnt of flanking DNA on each chrom. Rmap*/
#define DEFbeta  2.0        /*Beta parameter for effects Rqtl*/
#define DEFheritability 0.5 /*Heritability Rcross*/
#define SIG_LEVEL 3.84      /*Experimentwise sig. level Eqtl, Preplot*/
#define NUM_SIG 5           /*Number of background markers in CIM. Zmapqtl, JZmapqtl  (nbp)*/
#define WIN_SIZE 10.0       /*Window size in CIM. Zmapqtl, JZmapqtl */
#define DELTAX 2.0          /*Walking speed in cM.  Zmapqtl, JZmapqtl*/
#define MIN_DATA 25         /*Minimum sample size for typing markers. Qstats*/
#define HLINE     70        /*length of horizontal lines in output*/
/*  These are for the ECM algorithm*/
#define STOP_EM 0.000001
#define M_TIME 10000
#define MIN_DIST 0.0001  /*minimum distance from a marker for CIM analysis. */
#define MININTERVAL 0.0001
#define SIGNOFDEVIL 666.0  /* This just indicates that the ECM algorithm failed at a given site. */

/*These are for mapfunctions and what missing values are coded as */
#define MAPFUNCTION 1
#define MISS_VAL -1000000.0

/*Some general globals. */
#define MAXNAME 255
#define MAXLINE 511
#define SEED 11221960
#define BIG 2147483647.0   /*  = 2e31 - 1  */
#ifndef PI
#define PI 3.141592654
#endif
#define MAPDELTA 0.00001   /*Used in mapping functions for iterations  Utilities.c */

extern char time_buffer[MAXNAME+1];  /*String to keep the time */
extern char gbuffer[MAXLINE+1];      /*Reusable global buffer space  */
extern char gname[MAXNAME+1];        /*Reusable global name space    */
extern char gwarn[MAXNAME+1];        /*Reusable global warning space */

/* Some structures */

/*             QTL                  m3
c1 .......|.....K.....|.............|..........|....K....
   brdrs?              <-----s----->       
c2        |.....K.....|.............|..........|

brdrs are flanking DNA, that is DNA outside the first and last
markers on a chromosome.
*/

typedef struct MapofMarkers {
  int m;         /* number of chromosomes */
  int l;         /* average number of markers per chromosome */
  int maxl;      /* maximum number of markers on any chromosome */
  double sigl;   /* variance in number of markers per chromosome.  <= 0 => fixed */
  double s;      /* average distance between consecutive markers in cM */
  double sigs;   /* variance of s */
  int *mpc;      /* = ivector(1,m); number of markers for each chromosome.  = l forall iff sigl <= 0 */
  double **mrf;  /* = dvector(1,m,0,max(mpc)+1); If sigl <= 0, then l is max(mpc), and mpc isn't used. pointer to a matrix of recombination frequencies between markers i and i+1 */ 
  int ml;        /* total number of markers, = m*l iff sigl <= 0 */
  double brdrs;  /* Will there be chromosomal material outsite the flanking chromosomal markers? If this is negative, then there will be no borders.  */
  int *knum;     /* = ivector(1,traits) Number of QTLs on this map for each of the traits */
  int traits;    /* Number of traits to simulate */
  int otraits;   /* Number of other traits */
  int *otypes;   /* Number of classes for each of the other traits */
  char **tnames;  /* Names of the traits */
  char **onames;  /* Names of other traits */
  char **cnames; /* Names of the chromosomes */
  char **names;  /* Names of the markers */
  int  **ttable; /* Table to indicate where in names the marker name is */
  int  **types;  /* Indicate the type of marker:  
		    -1  =>  a-           a dom
		    0  =>  codominant
		    1  =>  A-           A dom           */
} markermap;

typedef struct QTLdata {
  int chrm;        /* this qtl is on chromosome chrm */
  int mrk;         /* it resides after this marker */
  double c1;       /* the recombination frequency between marker mrk and the qtl */
  double c2;       /* the recombination frequency between marker mrk+1 and the qtl */
  double a;        /* the additive effect of the qtl */
  double d;        /* the dominance effect of the qtl, 0.0 => no dominace */
  markermap *map;  /* pointer to the map of markers */
  int trait;       /* the trait that this is a qtl for. */
  double r2;       /*residuals from the estimation stage*/
  double tr2;      /*total residuals from the estimation stage*/
  double s;        /*test statistic for normality of residuals*/
}  aqtl;

typedef struct IndividualData {
  int **markers;  /* = gives the values of the markers */
  int **vqtls;    /* = imatrix(1,t,1,k) gives the alleles of the qtls */
  int t;          /* number of quantitative traits measured */
  double *y;      /* = dvector(1,t) gives the phenotypes for the t traits, 1 by default */
  char **oy;      /* other traits = cmatrix(1,ot,0,MAXNAME) */
  int *oyt;       /*  oyt = ivector(1,otraits)  */
  double *g;      /* gives the genotypic value over all QTLs */
  markermap *map; /* pointer to the map of markers */
  aqtl *qtls;     /* pointer to the array of qtl information */
  char *name;     /* name of this individual */
  char print_flag;/* should we print this individual? Or use it in analysis? */
  int bc;         /*which backcross in Design III.   not used yet.*/
}  individual;

typedef struct TotalGenome { /*Structure to hold a genome defined by the genetic linkage map*/
  int chrom;                 /* Chromosome of the marker interval */
  int markr;                 /* The marker that precedes the interval, can be 0 if borders are allowed */
  double dist;               /* This distance, from marker markr to markr+1 on chromosome chrm, is in M */
  double pos;                /* Position of marker from left telomere in Morgans*/
  int whichqtl;              /* Indicator of whether there is a qtl on the interval 1 if yes, 0 if no */
  double mxo;                /* xo position of maternal    xo = 0.0 if no crossover,   */
  double pxo;                /*             or paternal    xo > 0.0 the distance (M) from the marker to the xo */
  struct TotalGenome *prev;  /* Pointer to previous marker interval */
  struct TotalGenome *next;  /* Pointer to next marker interval */
}  genome;


typedef struct otraittype {   /*Structure to hold information about categorical traits and the categories*/
  char *name;                 /* Categorical (Otrait) trait name */
  int   which;                /* Which Otrait it is for */
  struct otraittype *prev;
  struct otraittype *next;
}  otnode;
  
typedef struct aline {   /*Structure to hold inbred lines that are simulated*/
  char *name;            /*Name of the inbred line*/
  individual *iptr;      /*pointer to the data for that line*/
  int nn;                /*sample size of this data set*/
  struct aline *prev;
  struct aline *next;
  int which;             /* */
}  thelines;

typedef struct Params {  /*Structure to hold parameters*/
  char *resource;           /*file that holds current parameter values */
  char *error;              /* log file*/
  char *map;                /* genetic linkage map*/
  char *mapin;              /* input file for Rmap*/
  char *qtl;                /* genetic model file*/
  char *eqtl;               /* file for estimates of genetic model*/
  char *qtlin;              /* Rqtl input file*/
  char *ifile;              /* data file */
  char *iinfile;            /* input for Rcross*/
  char *qstat;              /* output of Qstats*/
  char *lrfile;             /* output of LRmapqtl*/
  char *srfile;             /* output of SRmapqtl*/
  char *zfile;              /* output of Zmapqtl*/
  char *stem;               /* stem for filenames*/
  char *helpfile;           /* help file*/
  char *workdir;            /* working directory*/
  char *term;               /* terminal type (Preplot)*/
  char *tfile;              /* temporary file*/
  
  long seed;       /*random number seed*/
  char *thecross;  /*type of cross*/
  int chrom;       /*chromosomes in map*/
  int wchrom;      /*which chromosome to analyze*/
  int mark;        /*markers per chromosome*/
  double vmark;    /*st.dev. of markers/chrom*/
  double dist;     /*intermarker distance*/
  double vdist;    /*st.dev. of intermarker distances*/
  double tail;     /*flanking DNA*/
  int qnum;        /*number of QTL*/
  int dom;         /*dominance flag*/
  double beta;     /*parameter for effects simulation*/
  int traits;      /*number of traits*/
  int whichtrait;  /*which trait to analyze*/
  int reps;        /*number of repetitions in some simulations*/
  double Herit;    /*heritability*/
  double Environ;  /*envrinmental variance...overrides heritability*/
  int cross;       /*primary cross type */
  int crosst;      /*primary cross generations   */
  int tcross;      /*type of test cross*/
  int tcrosst;     /*generations of test cross*/
  int nn;           /*sample size*/
  int Model;        /*CIM analysis model*/
  double walk;      /*walking speed (cM) for CIM*/
  int Inter;        /*Interactive level flag*/
  int mapfunc;      /*mapping function*/
  int lodflag;      /*whether to convert to LOD scores*/
  double maxlr;     /*maximum LR.  Calculated by Eqtl and it ratchets.*/

  int gout;         /*output mode flag for Rmap*/
  int Rmode;        /*flag for map simulation type*/
  int verbosity;    /*flag to turn on/off verbose output*/
  int ihypo ;       /*Which hypothesis for Eqtl */
  int boot;         /*Prune flag for what to do*/
  double siglevel;  /*Experimentwise sig. level*/
  double size;      /*type I error probability*/
  double srf1;      /*alpha for adding markers in forward stepwise regression*/
  double srb1;      /*alpha for deleting markers in backward elimination regression*/
  int srm;          /*flag for which stepwise regression method to use*/

  int nbp;          /* # background markers in CIM*/
  double window;    /* window size in CIM*/
  int perms;        /*number of permutations in permutation test*/
  int boots;        /*number of bootstraps to do (actually, just a flag to do one)*/

  double null_sse;  /*Null hypothesis SSe*/
  double total_var; /*phenotypic sample variance of a trait*/
  double mapparam;  /*extra parameter needed by some mapping functions*/
} params;


#include "Utilities.h"
#if defined(MZFUNC)

typedef struct LinpakWorkspace {
  double **xx;       /* design matrix...dmatrix(1,p,1,n) */
  double **xsave;    /* save design matrix...dmatrix(1,p,1,n) */
  int ldx;           /* leading dimension of xx */
  int n;             /* columns of xx (sample size) */
  int p;             /* rows of xx  (background parameters) */
  int k;             /* should we make this the number of static rows?  mean plus otraits? */
  int t;                    /* number of traits in multitrait analysis */
  double **y;               /* trait matrix...dmatrix(0,t,1,n) */
  double *wy;              /* working trait vector...dvector(1,n) */
  double *qraux;     /* auxiliary vector...dvector(1,p) */
  double **rsd;             /* residual vector...dmatrix(0,t,1,n) */
  double **wrsd;            /* working residual vector...dmatrix(0,t,1,n) */
  double **bb;               /* regression coefficients vector...dmatrix(0,t,1,p) */
  double **jpvt;             /* pivot vector...dvector(0,t,1,p) */
  int **bp;          /* vector to indicate the chromosome and marker of background markers*/
  double *pp1;       /* a priori probability of QQ genotypes */
  double *pv;        /* a posteriori probability of QQ genotypes */
  double *pp2;       /* a priori probability of Qq genotypes */
  double *qv;        /* a posteriori probability of Qq genotypes */
  double **estimates;       /* parameter estimates for various hypotheses. 
			       = dmatrix(0,t,1,9)    */
  double **s2;       /* dmatrix(0,t,0,t) the variance-covariance matrix   */  
  double **s2i;      /* dmatrix(0,t,0,t) the inverse of the variance-covariance matrix   */  
  int *samplesize;   /* ivector(0,t) sample sizes for each trait.         */ 
  int *kpvt;         /*ivector(0,t) work space for invdet */
  double *work;      /*dvector(0,t) work space for invdet */
  int ipos;          /* the number of test positions */
  double *lratio;    /* dvector(1, ipos); Keep track of the LR of the sample */
  double *ts0;    /* dvector(0, t); null hypothesis Likelihood */
  int *pcnts;         /* ivector(1, ipos); Count number of shuffled sets with LR > SampleLR */
  
} linpakws;

linpakws *cd_linpakws(linpakws *ilnpk, int n, int p, int t, int cd);

#include "MZfunc.h"
#else

typedef struct LinpakWorkspace {
  double **xx;      /* design matrix...dmatrix(1,p,1,n) */
  double **xsave;   /* save design matrix...dmatrix(1,p,1,n) */
  int ldx;          /* leading dimension of xx */
  int n;            /* columns of xx */
  int p;            /* rows of xx */
  int k;
  double *y;     /* trait vector...dvector(1,n) */
  double *wy;
  double *qraux; /* auxiliary vector...dvector(1,p) */
  double *rsd;   /* residual vector...dvector(1,n) */
  double *wrsd;
  double *bb;    /* regression coefficients vector...dvector(1,p) */
  int *jpvt;  /* pivot vector...dvector(1,p) */
  int **bp;
  double *pp1;
  double *pv;
  double *pp2;
  double *qv;
  double *estimates;
  double *thestats;
} linpakws;

linpakws *cd_linpakws(linpakws *ilnpk, int n, int p, int cd);

#endif


#endif








