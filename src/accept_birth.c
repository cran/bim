/****************************************************
  File:       accept_birth.c
  Written by: Patrick Gaffney
  Date:       November 11, 2000
  Version:    0.4

  Purpose:

    Function ACCEPT_BIRTH computes the acceptance
    probability for the birth of a new QTL.
****************************************************/

#include "revjump.h"
#include "ranlib.h"

void checkResid(int nn, int nQtl, double mu, double* y, 
		QTL_INFO** allQtls, double* resid,int gmiss)
{
  double sum;
  int i,j;
  int* geno;
  QTL_INFO* aqtl;

  for(i=1;i<=nn;i++){
    sum = mu;
    for(j=1; j<=nQtl; j++) 
      {
	aqtl = allQtls[j];
	geno = igenotype(aqtl);
	sum += EFFECT_VALUE(geno[i], aqtl->a);
      }
    sum = y[i] - sum;
    if (!EQUALS(resid[i], sum)) 
      {
  for(i=1;i<=nn;i++){
    sum = mu;
    for(j=1; j<=nQtl; j++) 
      {
	aqtl = allQtls[j];
	geno = igenotype(aqtl);
	Rprintf("%d(%lf) ",geno[i],  EFFECT_VALUE(geno[i], aqtl->a));
	sum += EFFECT_VALUE(geno[i], aqtl->a);
      }
    Rprintf("\n");
  }
	exit(1);
	Rprintf("==>%d .. %f %f",i,sum, resid[i]);
      }
  }
}



int select_move(int nQtl, double* bp, double* dp)
     /*---------------------------------------------------------------
 * Description
 *    Selects the move type (1=>birth, 2=>death, 3=>none)
 *
 * chrom           pointer to chromosome where death/birth is    
 *                 occuring                                      
 * qtl_prior_mean  prior mean for the number of qtl (for POISSON!)
 * cval            Multiplier for birth and death probabilities.
 *                 Must be chosen so that the product cval*
 *                 (birth_prob+death_prob)<1. Recommended cval=0.4
 *                 If cval is in error, the maximum value is chosen
 *                                                               
 *--------------------------------------------------------------*/								    
{
  double uni_ran;
  int move_type;
  


  uni_ran = genunf(0.0, 1.0);
  if( uni_ran < bp[nQtl])
    move_type = C_BIRTH;
  else if( uni_ran < (bp[nQtl] + dp[nQtl]) )
    move_type = C_DEATH;
  else
    move_type=C_UPDATE;


  return move_type;
}





double get_log_proposal_ratio(int nQtl, double* bp, double* dp, double* priorRatio)
{
  /* nQtl ... number of QTL in reduced model */
  return log( (dp[nQtl+1]/bp[nQtl]) * priorRatio[nQtl+1]);
}



double get_log_position_ratio(int revjump, CHROMOSOME* chrom, DATA* myData)
{
  if (revjump & MCMC_FLAG_RANDOM_CHROM)
    return log(chrom->chromLen * myData->nChrom / myData->totalChromLen);
  else
    return 0;   /* we sample the position from the prior */
}

void calcResid2(int nn, int nQtl, double* y, double* modified_parms,
		QTL_INFO** qtls, double* newResid)
{
  int i,j,idx;
  QTL_INFO* aqtl;
  double a[QTL_ADD_DOM];
  int* geno;

  for(i=1;i<=nn;i++) newResid[i] = y[i] - modified_parms[1];

  for(j=1, idx=2; j<=nQtl; j++) 
    {
      aqtl = qtls[j];
      geno = igenotype(aqtl);

      /* find which terms are applicable for the QTL */
      a[QTL_ADD] = (aqtl->flag & QTL_ADD)? modified_parms[idx++]: 0;
      a[QTL_DOM] = (aqtl->flag & QTL_DOM)? modified_parms[idx++]: 0; 

      for (i=1;i<=nn;i++) 
	newResid[i] -= EFFECT_VALUE(geno[i], a);
    }
}







