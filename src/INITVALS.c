/**************************************************************
  File:       initgeno.c
  Written by: Patrick Gaffney
  Date:       November 11, 2000
  Version:    0.4

  Purpose:
  -------
    These functions find the index of missing markers and traits.

**************************************************************/

#include "revjump.h"
#include "ranlib.h"

/********************************************
         Function INITQTL_GENO

  initializes qtl genotypes given their 
  location and flanking marker genotypes
  for the given chromosome
********************************************/


void initQtl(int nn, QTL_INFO* aqtl, 
	     params* theparams, int bc, int gmiss, int revjump)
{
  /* here we initialize the QTL on the basis of its flanking markers 
     only (and don't use trait values as in propose_qtl_geno() in
	file update_qtl.c) 
*/
	int i;
	int* geno;
	mygenome* qptr = aqtl->qptr;

    genProbs(nn, theparams, bc, gmiss, aqtl, revjump);   
    aqtl->prevDist = qptr->prev->dist;
    aqtl->nextDist = qptr->dist;
	
	
	geno = igenotype(aqtl);

	for (i=1; i<=nn; i++)
	   GenGenotype(gmiss, aqtl->log_prob[i], &geno[i]);
}




void GenGenotype (int gmiss, double* log_pr, int* geno)				 
     /* assumes pr indexed by -1,0 or 1 */
     /* gmiss is set only for crosses with only 2-genotype progeny, to
   indicate the missing genotype (-1,0 or 1) */
{
  double u;
  
  u = genunf(0,1);	
  
  if (log(u) < log_pr[1] && gmiss !=1) *geno = 1;
  else if (log(1-u) <= log_pr[-1] && gmiss != -1) *geno = -1;
  else *geno = 0;
}




int* igenotype(QTL_INFO* aqtl)
/* the first element (zero index) of array contains the first genotype value */
{
	return aqtl->qptr->genotype;
}


int** p_igenotype(QTL_INFO* aqtl)
/* the first element (zero index) of array contains the first genotype value */
{
	return &aqtl->qptr->genotype;
}


int* genotype1(QTL_INFO* aqtl)
/* the first element (zero index) of array contains the first genotype value */
{
	return (&aqtl->qptr->genotype)[1];
}



void dgenotype(QTL_INFO* aqtl, int nn, double* dz)
/* the first element (zero index) of array dz contains the first genotype value */
{
	int i;
	int* geno = igenotype(aqtl);

	for (i=0; i<nn; i++)
		dz[i] = (double)(geno[i]);
}



int qtlnum(QTL_INFO* aqtl)
{
	return -aqtl->qptr->markr;
}



double IMdist(QTL_INFO* aqtl)
{
	CHROMOSOME* chrom = aqtl->chrom;
	int mark = binSearch(chrom->nMark, chrom->mark_pos, aqtl->qptr->pos);
	return chrom->mark_pos[mark+1] - chrom->mark_pos[mark];
}


/*************************************************
             FUNCTION HALDANE
*************************************************/

double haldane( double dist_val ){
/* expects dist_val in Morgans */
  double recval;
  recval = 0.5 * ( 1 - exp( -2.0 * dist_val ) );
  return recval;
}



void random_perm_QTL(int nQtl, QTL_INFO** qtls, QTL_INFO** buff, int* perm_num)
/* randomly permute the QTLs (for updating) */
{
	int i,j, val;

	for (i=1; i<=nQtl; i++) perm_num[i] = i;
    for (i=nQtl; i>=1; i--) 
      {
        j = ignuin(1,i);
	val = perm_num[j];

	SwapInt(&perm_num[i], &perm_num[j]);  /* now move entry containing val to end */
	buff[i] = qtls[val];
	buff[i]->qptr->markr = -i;            /* make sure marker numbers consistant */
	/* though that now doesn't mean much   */
      }
    for (i=1; i<=nQtl; i++) 
      qtls[i] = buff[i];                       /* restore our list */
}



void printX(int nn, int nQtl, QTL_INFO** qtls)
{
  FILE* f;
  QTL_INFO* aqtl;
  int i,j;

  f = fopen("x.txt","w");
  for (i=1; i<=nn; i++)
    {
      fprintf(f,"1  ");
      for (j=1; j<=nQtl; j++)
	{
	  aqtl = qtls[j];
	  if (aqtl->flag & QTL_ADD) fprintf(f,"%d  ",aqtl->qptr->genotype[i]);
	  if (aqtl->flag & QTL_DOM) 
	    fprintf(f,"%5.2f  ",DOM_GENO_VALUE(aqtl->qptr->genotype[i]));	  
	}
      fprintf(f,"\n");
    }
  fclose(f);
}

void printXtX(int nterm, double** XtX)
{
  int i,j;
  FILE* f;
	  
  f = fopen("x.xtx","w");
	  
  for (i=1; i<=nterm; i++)
    {
      for (j=1; j<=nterm; j++)
	fprintf(f,"%f ",XtX[i][j]);
      fprintf(f,"\n");
    }
  fclose(f);
}



double gammln2(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}
