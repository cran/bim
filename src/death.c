/****************************************************
  File:       death.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:

    Function DEATH performs the "death step" for
    removing an existing QTL.
****************************************************/

#include "revjump.h"
#include "ranlib.h"




int death(DATA* myData, PRIORS* priors, 
		   params* theparams, WORK* myWork, MCMC_PARAM* myMCMC)
{
  double log_accept_death_prob, log_effect_ratio, 
	     log_proposal, log_position;
  double uni_ran;

  /* get elements */
  int nQtl = myData->nQtl;
  int nn = myData->nn;
  double mu = myData->mu;
  double sigmasq = myData->sigmasq;
  QTL_INFO** all_qtls = myData->myQtls;
  QTL_INFO* deathQtl;
  CHROMOSOME* chrom;
  double bCb, d_invD_d, log_det_Chol, log_det_invD;

  /* work structures we'll use */
  double* mod_effect = myWork->mod_effect;    /* new grand mean and effects */
  double** XtX = myWork->XtX;                 /* XtX matrix */
  double* XtY = myWork->XtY;                 /* XtX vector */
  double** chol = myWork->chol;               /* new Cholesky decomposition */
  double* p = myWork->p;                      /* new diagonal for the Cholesky */
  double* pmean = myWork->pmean;              /* mean of grand mean/effects */
  double* pvar = myWork->pvar;                /* variance of grand mean/effects */
  double* u = myWork->u;                      /* iid normal values */
  double* weight = myWork->weight;            /* new weights */
  int nterm = myData->na[ADD] + myData->na[DOM] + 1;
    int changeQtl;



  changeQtl = ignuin(1,nQtl);
  if (changeQtl != nQtl)
     moveQtlToEndOfXtX(nQtl, all_qtls, changeQtl, XtX, XtY, nterm);					 

  deathQtl = all_qtls[nQtl];   /* recall we've randomized order */
  chrom = deathQtl->chrom;

  log_position = get_log_position_ratio(myMCMC->revjump, deathQtl->chrom, myData);

  log_effect_ratio = get_effect(nn, nQtl, myData->y, mu, sigmasq,
	                            deathQtl, NULL, all_qtls, myMCMC->revjump,
			                    priors, myData->na, mod_effect,
  			                    XtX, XtY, chol, p, pmean, pvar, u, weight,
								&bCb, &d_invD_d, &log_det_Chol, &log_det_invD, myWork);


                     
  log_proposal = get_log_proposal_ratio(nQtl-1, myData->bp, myData->dp, myData->prior_ratio);

  log_accept_death_prob = -log_proposal - log_position + log_effect_ratio;


/* DECIDE WHETHER OR NOT TO ACCEPT NEW QTL */

/**************************************************
proceed with accept/reject death only if jacobian
is not 0. if jacobian is 0, reject death.
**************************************************/


    uni_ran = genunf(0.0,1.0);
    uni_ran = log(uni_ran);

    if(uni_ran < log_accept_death_prob) 
	{
      dropQtl(nn, &myData->nQtl, myData->y, all_qtls, mod_effect, weight, 
		      &myData->mu, myWork->resid, myData->na);

#ifndef MCMC_OMIT_DEBUG
	 checkIntegrity(myData->nQtl, deathQtl->chrom);
#endif

	 setCholParams(myWork, bCb, d_invD_d, log_det_Chol, log_det_invD);
	  return 1;
	}
	return 0;
}




void dropQtl(int nn, int* nQtl, double* y, QTL_INFO** all_qtls, 
	         double *mod_effect, double* w, double* mu, double* resid, int* na)
{
  /*  recall  mu = effect[0]  and modified_mu = modified_parms[0]
      also modified_parms[*nqtl+1] contains the effect of the new qtl */

  QTL_INFO* dropQtl = all_qtls[*nQtl];
   
  /* remove its record from genome list */
  removeQtl(dropQtl->qptr);
  if (*nQtl == 0)
	  printf ("cannot delete non-existant QTL\n");
  else
     (*nQtl)--;

 
  /* set values */
  setEffect(nn, *nQtl, y, mod_effect, all_qtls, mu,
	        w, resid, na);
  dropQtl->chrom->nQtl--;
}




/***************************************
      Function GET_DEATH_INTERVAL
    randomly draws a new QTL interval 
***************************************/

int get_death_interval( int nQtl )
{
  return (int)((genunf(0.0, 1.0) * nQtl) + 1);    /* recall qtls are numbered 1,2,... */
}





int removeQtlFromList(int nQtl, QTL_INFO* aqtl, QTL_INFO** qtls)
{
	int i,j;

	for (i=1; i<= nQtl; i++)
	{
		if (qtls[i] == aqtl)
		{
			for (j=i+1; j<=nQtl; j++) 
				qtls[j-1] = qtls[j];
			break;
		}
	}
	return i;
}





