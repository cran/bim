#ifndef QTL_GUARD
#define QTL_GUARD
/**************************************************************
  File:         qtl.h
  Written by:   Patrick Gaffney
  Modified by:  Amy Jin
  Date:         November 11, 2000/Aug 13,2002
  Version:      0.4

  Purpose:
  -------

    Header file describing qtl records in MCMC program.  

**************************************************************/

#include "mygenome.h"
#include "chrom.h"

/* codes which tell which effects a QTL has 
  (C_NONE is generally when there are
   only interactive terms for the QTL 
*/
#define QTL_NONE    0
#define QTL_ADD     1
#define QTL_DOM     2
#define QTL_ADD_DOM 3   /* QTL_ADD | QTL_DOM */


#define EFFECT_VALUE(x, a) ((x)==0? (a)[DOM]: (a)[ADD]*(x)) 
#define DOM_GENO_VALUE(x) ((x)==0? 1:0)
#define DOM_VALUE(x,dom) ((x)==0? (dom):0)
#define ADD_VALUE(x,add) (x)*(add)
#define IGENO_VALUE(x,type) ((type)==ADD ? (x): ((x)==0? 1: 0))
#define GENO_VALUE(x,type) ((type)==ADD ? (double)(x): ((x)==0? 1.0: 0.0))


/*
#define EFFECT_VALUE(x, a) ((x)==0? (a)[DOM]*0.5: (a)[ADD]*(x)-(a)[DOM]*0.5) 
#define DOM_GENO_VALUE(x) ((x)==0? 0.5:-0.5)
#define DOM_VALUE(x,dom) ((x)==0? (dom)*0.5:-(dom)*0.5)
#define ADD_VALUE(x,add) (x)*(add)
#define GENO_VALUE(x,type) ((type)==ADD ? (double)(x): ((x)==0? 0.5: -0.5))
*/


/* alternative coding 
========================================================================
#define EFFECT_VALUE(x, a) ((x)==0? a[QTL_DOM]: a[QTL_ADD]*x) 
#define GENO_VALUE(x,type) (type==QTL_ADD ? (double)x: (x==0? 1.0: 0.0))
========================================================================
*/


/* typedef struct chromosomeData CHROMOSOME; */
struct MCQtlData
{
  CHROMOSOME* chrom;   


  mygenome* qptr;      /* pointer to genome record, 
                        NOTE: this contains number of QTL
			      (order it appears) on chromosome
                      */

  
  int flag;           /* QTL_NONE, QTL_DOM, QTL_ADD, QTL_ADD_DOM */
  int nParam;         /* number of parameters (determinable from flag), but
                         handy to have around */
  int nInt;           /* number of interaction terms */
  double* a;          /* effects ... a[QTL_ADD] is additive,
                                     a[QTL_DOM] is dominance */
  double* w;          /* weighting for prior variance (1/p)
					             ... w[QTL_ADD] is additive,
                                     w[QTL_DOM] is dominance */


  /* table of conditional probabilities of QTL given flanking markers */
  double** transProbL;  /* transition probabilities to left marker */
  double** transProbR;  /* transition probabilities to right marker */
  double*** log_condProb;
  double** log_prob;    /* a prob table table for genotype given marker for each
                           individual. Possibly useful for missing data 
                           In reality, this is an array of pointers, either
                           to elements (3x1 vectors) of log_condProb (when flanking
						   markers are fully informative) or to the corresponding
                           element of missing_prob otherwise.  */
  int probValid;        /* a flag indicating the 'prob' is still valid... used
						   in updates on one QTL */


  double** missing_prob; /* */

  double prevDist;   /* used for determining if transProbL should be altered */
  double nextDist;   /* used for determining if transProbR should be altered */ 
  
  
};



#endif






