/**************************************************************
  File:         read_data.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:
  -------
    getdata of type VOID.This reads 
    trait, marker, distance and initial values from different 
    files. 



**************************************************************/

#include "revjump.h"
#include "ranlib.h"

/*************************************************
             FUNCTION GETDATA
*************************************************/

void getdata( MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, 
              CHROMOSOME* chromInfo, WORK* myWork,
	      int* revjump, double* burn_in, double* preburn_in, 
	      int* niter, int* nby, double* mu, double* sigmasq,
	      double* meanmu, double* sdmu, double* siga1, double* siga2,
	      double* meanadd, double *sdadd,double* meandom, double*sddom,
	      double* qtlmean)
{

  double temp;
  double prior_sd[3];
  int type;


  myMCMC->nby = *nby;
  myMCMC->niter=*niter;
  myMCMC->revjump=*revjump;
  myData->mu=*mu;
  myData->nQtl=0;
  myData->sigmasq=*sigmasq;
  priors->mean[MU]=*meanmu; 
  prior_sd[MU]=*sdmu;
  priors->sig_a1=*siga1; 
  priors->sig_a2=*siga2;
  priors->mean[ADD]=*meanadd;
  prior_sd[ADD]=*sdadd;
  priors->mean[DOM]=*meandom;
  prior_sd[DOM]=*sddom; 
  
  

   /*------------------------------------------------------------------ 
   * determine how many parameters to add each time in birth/death step
   * if revjump field is set for RJ_MCMC, we only have birth/death 
   * of one parameter at a time.  So addParam = 1.  If gmiss is not -2,
   * then we have only 2 genotypes, so addParam = 1.  Otherwise, we
   * will favor adding an additive and dominance parameter (addParam=2)
   *------------------------------------------------------------------ 
   */	  
  myMCMC->burnIn = (*burn_in >= 0 && *burn_in < 1)? *burn_in: -1;
  myMCMC->preBurnIn = (*preburn_in >= 0 && *preburn_in < 1)? *preburn_in: -1;


  if (myData->gmiss == -2 && 
	  (myMCMC->revjump & MCMC_SELECT_DOM_FIELD) != SELECT_RJMCMC &&
	  ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_GET_APPROP))
	   myMCMC->addParam = QTL_ADD_DOM;
  else
  {
	  if ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_GET_DOM_ONLY)
	     myMCMC->addParam = QTL_DOM;
	  else if ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_RJMCMC && myData->gmiss == -2)
	     myMCMC->addParam = QTL_NONE;
	  else
	     myMCMC->addParam = QTL_ADD;
  }

  
/*************************************************************************
   the parameter revjump is binary. revjump=1 means do 
   reversible jump mcmc computations. if revjump=0, no reversible
   jump mcmc, do only regular mcmc as in satagopan et al (1996).
*************************************************************************/

 
  /* initialize our data structures */
  if (priors->sig_a2 < 0) priors->sig_a2 = -priors->sig_a2 * myData->y_var;
  
  
  /* qtl prior mean is required only if doing reversible jump mcmc
     computations as the following "if" statement indicates */
  
  if( (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD) >= 1) {
    priors->qtl_mean= *qtlmean;
   }
  else{
    priors->qtl_mean = myData->nQtl;
  }
  
  
  /* check for defaults (-1 in standard deviation field => 
     use the mean field (0<mean<1 typically) times the data
     variance to give the variance for the prior variance */
  if (myData->sigmasq<0) 
    {
      myData->sigmasq = myData->mu * myData->y_var;
      myData->mu = myData->ybar;
    }
  
  priors->sampleVar[MU] = 0;
  
  if (prior_sd[MU] < 0) 
    {
    priors->var[MU] = priors->mean[MU] * myData->y_var;
    priors->mean[MU] = myData->ybar;
  }
  else
    priors->var[MU] = prior_sd[MU]*prior_sd[MU];
  
  priors->log_var[MU] = log(priors->var[MU]);
  
  
  
  for (type = ADD; type <= DOM; type++)
    {
      priors->sampleVar[type] = 0;
      priors->alpha[type]= priors->beta[type]=0.0;
    
    if (LE(prior_sd[type],0)) 
      {
	if (!LE(priors->mean[type],0)) 
	  priors->var[type] = priors->mean[type] * myData->y_var;
	else
	  {
	    priors->sampleVar[type] = 1;
	    
	    if (EQUALS(prior_sd[type],0) || !LT(priors->mean[type], 0))
	      {
		priors->alpha[type] = 2;
		priors->beta[type] = 10;   /* this gives Beta(1,5) as a prior with mean 0.25 */
	      }
	    else 
	      {
		priors->alpha[type] = -priors->mean[type];
		priors->beta[type] = -prior_sd[type];
	      }
	    temp= genbet(priors->alpha[type], priors->beta[type]);
	    priors->var[type] = temp * 2 * myData->y_var;	
	  }
	priors->mean[type] = 0;
      }
    else
      priors->var[type] = prior_sd[type]*prior_sd[type];
    
    priors->log_var[type] = log(priors->var[type]);
    }
}




