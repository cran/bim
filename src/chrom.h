#ifndef CHROM_GUARD
#define CHROM_GUARD
/**************************************************************
  File:         chrom.h
  Written by:   Patrick Gaffney
  Modified by:  Amy Jin
  Date:         November 11, 2000/Aug 13,2002
  Version:      0.4

  Purpose:
  -------

    Header file describing chromosome records in MCMC program. 



**************************************************************/


#include "mygenome.h"


typedef struct chromosomeData CHROMOSOME; 
typedef struct MCQtlData QTL_INFO;


struct chromosomeData
{
   int num;
   int nQtl;
   int nMark;             /* number of markers in the chromosome */
   double chromLen;       /* length of the chromosome (in Morgans) */

   mygenome** mark_genome;  /* array of pointers to the genome entries for markers */
   double* mark_pos;
  
   QTL_INFO** qtls; 

   double** log_prob;

};

#endif




