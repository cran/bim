/************************************************************
    File:       mcmc.c
    Written by: Patrick Gaffney
    Date:       November 11, 2000
    Version:    0.4

    Purpose:

      Function MCMC gets the states of the required 
      Markov chain from its equilibrium distribution.
************************************************************/
#define PRE_BURN_IN .05   /* accelerated percent of burn-in     */
#define BURN_IN .05      /* burn-in is BURN_IN percent of total*/
#define ELEMENTS_PER_QTL_OUTPUT 13    /* number of elements in a qtl output "line" */


#include "revjump.h"
#include "ranlib.h"
#include <assert.h>


void mcmc(MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors,
          CHROMOSOME* chrInfo,
          WORK* myWork, double* result) 
{
  int move_type, num_qtl[MAXQTL+1];
  int output_idx, allocRec;

  int iter, i;
  int nn = myData->nn;
  double** nmove = dmatrix(1,5,1,2);
  QTL_INFO** all_qtls = myData->myQtls;
  int countedIter=0;
  int always_update = 1;                   /* for pre-burn-in */  
  /* a double negative, forces update of peramaters
     regardless of whether birth/death ocurred */
  int always_update_long = (myMCMC->revjump & MCMC_OPTIMIZE_BOOST);
  int rj_mcmc = (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD);
  QTL_INFO* aqtl;
  int k;

  /*write_res = fopen("test", "w");
  fprintf(write_res, 
  "niter nqtl iqtl chrom LOD mu sigmasq addvar domvar add dom locus esth");
  if (!rj_mcmc) fprintf(write_res, " HM SHM\n");
  else fprintf(write_res, "\n");
  */

  /* index for writing to result */
  output_idx=0;
  allocRec=0;

  /* adjust to get posterior */
  priors->sig_a1 += nn / 2.0;
  initVars(nn, nmove, num_qtl, myMCMC, priors->qtl_mean, 1.0);
  setupBurnIn (myMCMC, myData);
	               
  
  
  setupBirthDeathProbs(SELECT_NQTL_UNIFORM, priors->qtl_mean, 
		       0.9, &myData->bp, &myData->dp, &myData->prior_ratio);




  for( iter=(int)-myMCMC->burnIn; iter < myMCMC->niter; iter++)
    {

      /* restore the user's choice */
      if (iter == (int)(-myMCMC->burnIn + myMCMC->preBurnIn))
	{
	  if (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD)
	    setupBirthDeathProbs(myMCMC->revjump, priors->qtl_mean, 
				 myMCMC->cval, &myData->bp, &myData->dp, 
				 &myData->prior_ratio);
	  always_update = (myMCMC->revjump & MCMC_SELECT_100PCNT_UPDATE);
	}



     if (iter % 5000 == 0)  /* perform integrity checks */
	{
	  Rprintf ("%d.",iter/1000);
	  checkResid(myData->nn, myData->nQtl, myData->mu, myData->y, 
		     myData->myQtls, myWork->resid, myData->gmiss);
	  for (i=1; i<= myData->nChrom; i++)
	    checkIntegrity(myData->nQtl, &chrInfo[i]);
	} 

      /****************************************
       *determine whether birth or death step
       *or whether to update parameters only
      ****************************************/

      /*********************************************************************
      **********************************************************************
      BIRTH AND DEATH STEPS ARE REQUIRED ONLY WHEN DOING REVERSIBLE
      JUMP MCMC COMPUTATIONS, AS THE FOLLOWING IF STATEMENT INDICATES.

      FOR REVERSIBLE JUMP MCMC, CHOOSE MOVE TYPE (BIRTH, DEATH OR
      STAY) AND THEN PROCEED. 
      **********************************************************************
      *********************************************************************/


      if (rj_mcmc) 
	{
	  move_type = select_move( myData->nQtl, myData->bp, myData->dp);
	  if (move_type != C_UPDATE)
	    birth_death(move_type, myData, myMCMC, myWork, priors, chrInfo, nmove);
	}
      else
	move_type = C_UPDATE;           /* always update for fixed MCMC */


      if (always_update) move_type = C_UPDATE;        /* force 100% updates */
		

	

      if (move_type == C_UPDATE)
	fixed_locus_update(nn, myData->nQtl, all_qtls, myMCMC, myData, 
			   priors, myWork, myData->theparams, nmove);
      else if (myData->nQtl >=1 && always_update_long)
	{
	  nmove[C_EXCHANGE][1]++;
          
	  aqtl = long_range_update(nn, myData->nQtl, all_qtls,
				   myData->chrom_pos, myData->totalChromLen, 
				   myData->nChrom, myData->myChroms,
				   priors, myData, myWork, myMCMC);	
	  if (aqtl)
	    {
	      nmove[C_EXCHANGE][2]++;
	      setValidFlag(aqtl, myMCMC->revjump);         
	    }
		  
	}



      /* record frequency of the number of QTL after burn-in (every entry) */
      if (iter >= 0) num_qtl[myData->nQtl]++;

      /* output results and possibly perform diagnostics */
      if(iter % myMCMC->nby == 0)
	{
	  countedIter++;
	  outputResults(nn, iter, myData, all_qtls, myWork->resid, 
			myWork->output_qtls, priors, myMCMC, result, &output_idx);
	}


    }  /* end iter */

  /* close results file 
  fflush(write_res);
  fclose(write_res);*/
}



void birth_death(int move_type, DATA* myData, MCMC_PARAM* myMCMC, WORK* myWork, PRIORS* priors, 
		 CHROMOSOME* chrInfo, double** nmove)
{
  static char* moves[] = {"-","B","D","U"};


  nmove[move_type][1] += 1;

  if (move_type==C_BIRTH)
    {
      if (birth(myData, chrInfo, priors, myData->theparams, myWork, myMCMC))
	nmove[C_BIRTH][2] += 1;
    }
  else if (move_type==C_DEATH) 
    {
      if (death(myData, priors, myData->theparams, myWork, myMCMC))
	nmove[C_DEATH][2] += 1;
    }

  if (myData->nQtl >= 1 && 
      ((myMCMC->revjump & MCMC_SELECT_DOM_FIELD) == SELECT_RJMCMC))
    {
      nmove[C_SWAP][1]+=1;
      if (swap_add_dom(myData->nn, myData->nQtl, 
		       myData->myQtls, priors, myData, myWork, myMCMC))
	nmove[C_SWAP][2] += 1;
    }
				   

}




void fixed_locus_update(int nn, int nqtl, QTL_INFO** all_qtls,MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, WORK* myWork, params* theparams, double** nmove)
{
  int i; 
  QTL_INFO* aqtl;
  double sum;
  int simpleGibbs = myMCMC->revjump & MCMC_SELECT_SIMPLE_GIBBS;

	  

  if(nqtl > 1) 
    random_perm_QTL(nqtl, myData->myQtls, myWork->output_qtls, myWork->perm_num);

  /* we must update the sigmasq and effect_prior_var BEFORE update_effects 
     to retain the validity of the Cholesky decomposition 
  */
  /*  myData->sigmasq = MH_update_sigmasq(nn, myWork->resid,
      myData->sigmasq, 
      myData->y_var/20.0, 1.5 * myData->y_var);
  */
  myData->sigmasq = Gibbs_update_sigmasq(nn, myWork->resid, priors->sig_a1, priors->sig_a2);

  if (priors->sampleVar[ADD] || priors->sampleVar[DOM])
    {
      if (myMCMC->addParam == QTL_NONE || myMCMC->addParam == QTL_ADD_DOM)
	{
	  if (priors->sampleVar[ADD]) 
	    update_effect_prior_var(ADD, nqtl, all_qtls, priors, myData->y_var);
	  if (priors->sampleVar[DOM]) 
	    update_effect_prior_var(DOM, nqtl, all_qtls, priors, myData->y_var);
	}
      else if (priors->sampleVar[myMCMC->addParam])
	update_effect_prior_var(myMCMC->addParam, nqtl, all_qtls, priors, myData->y_var);
    }


  for (i=1, sum=0; i<=nqtl; i++) 
    {
      aqtl = all_qtls[i];	
      sum += update_lambda_qtl(aqtl, all_qtls, myData, myWork, myMCMC->revjump);


     if (simpleGibbs)
	update_effect(nn,nqtl, aqtl, all_qtls, myMCMC, myData, priors, myWork);	
      setValidFlag(aqtl, myMCMC->revjump);         
	                        /* invalidate the cache of log_prob of QTL given flanking
				   markers, for safety if warranted (i,e. more than one
				   QTL on this chromosome, or if feature disabled) */

    }


  if (!simpleGibbs)
    {
      update_effects(nn, nqtl, myMCMC->revjump, 
		     myData->y, &myData->mu, myData->sigmasq, myData->na,
		     all_qtls, priors, myWork);

    }
  else
    myData->mu = update_mu(nn, myData->mu, myData->sigmasq, myWork->resid,
			   priors->mean[MU], priors->var[MU]);



  /* now propose long-range updates (if appropriate) */
  if (nqtl >= 1 && 
      !(myMCMC->revjump & MCMC_SELECT_MOVE_NONE)) 
    {	
      nmove[C_EXCHANGE][1]++;
      aqtl = long_range_update(nn, nqtl, all_qtls, 
			       myData->chrom_pos, myData->totalChromLen, 
			       myData->nChrom, myData->myChroms,
			       priors, myData, myWork, myMCMC);					  
      if (aqtl)
	{
	  nmove[C_EXCHANGE][2]++;
	  setValidFlag(aqtl, myMCMC->revjump);         

	}


    }

		  
  /* keep tally */
  if (nqtl >= 1) 
    { 
      nmove[C_UPDATE][1]++;
      nmove[C_UPDATE][2] += sum/myData->nQtl;
    }

}






void update_effect(int nn, int nqtl, QTL_INFO* aqtl, QTL_INFO** all_qtls,
		   MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, WORK* myWork)
{
  if (aqtl->flag & QTL_ADD)
    aqtl->a[QTL_ADD] = update_add_effect(nn, myData->nQtl, myData->sigmasq, aqtl, 
					 myWork->resid, priors->mean[QTL_ADD], 
					 aqtl->w[QTL_ADD] * priors->var[QTL_ADD]);

  if (aqtl->flag & QTL_DOM)
    aqtl->a[QTL_DOM] = update_dom_effect(nn, myData->nQtl, myData->sigmasq, aqtl, 
					 myWork->resid, priors->mean[QTL_DOM], 
					 aqtl->w[QTL_DOM] * priors->var[QTL_DOM]);
}


int numcmp(double *v1, double *v2)
{
  if(*v1<*v2)
    return -1;
  else if(*v1>*v2)
    return 1;
  else
    return 0;
}


int binSearch(int nval, double* vals, double searchVal)
     /* assumptions:
      searchVal lies between vals[1] and vals[nval]

    inputs:
	    nvals    :  number of intervals
		vals     :  array containing nval interval bounds
		searchVal:  value whose interval we're seeking
    returns:
	  the interval in which searchVal lies, i.e. 
	      i : vals[i]<=searchVal<vals[i+1]
*/
{
  int low, high, mid;
  low = 1;
  high = nval;
  do {
    mid = (low+high)/2;
    if (searchVal < vals[mid]) {
      if (high==mid) 
	break; 
      high=mid;
    }
    else {
      if (low==mid) break; 
      low=mid;
    }
  } while (1);

  return low;
}


double lodnull(int nn, double y_var)
{
  return -(nn/2)*log(y_var) - nn/2;
}



double get_lod(int nn, double sigmasq, double ybar, double y_var,
	       double* resid)
{
  double lod, null_lod, work, sum;
  static double logten = 2.302585;
  int i;

  sum = 0.0;
  for (i=1; i<= nn; i++) {
    work = resid[i];
    work *= work;
    sum += work;
  }

  lod = -(nn/2)* log(sigmasq) - 0.5*sum/sigmasq;
  null_lod = lodnull(nn, y_var);

  return (lod - null_lod)/logten;
}


void normalProb(int nn, double* resid, double sigmaSq, double* diag_fY,
		double* maxResid, double* pcnt_gt_99)
{
  /* increment the f(y) diagnostic by */

  static double logTwoPi = 1.837877066;
  double sqrtSigmaSq = sqrt(sigmaSq);
  int i;
  double z;
  double zlim;
		
  if (40*nn < 4000)	zlim = 1.412979 * pow(40.0*nn,0.1122);   
	                        /* approx to find z such that P(z<-zlim or z>zlim) = p = 1/(20*nn)
	                           reasonable for 0.01 < p < 0.0005 */
  else zlim = 2.221161 * pow(40.0*nn, 0.05646);
	                        /* approx to find z such that P(z<-zlim or z>zlim) = p = 1/(20*nn)
	                           reasonable for 0.0005 < p < 0.000001 */
  if (zlim < 2) zlim=2;

  for (i=1; i<=nn; i++)
    {
      z = resid[i]/sqrtSigmaSq;
      if (fabs(z) >  fabs(maxResid[i]))  maxResid[i] = z;

		
      if (fabs(z) > zlim) pcnt_gt_99[i] += 1.0;

      diag_fY[i] += exp(-logTwoPi -resid[i]*resid[i]/(2*sigmaSq))/sqrtSigmaSq;
    }
}







/* results are output as follows
     for each QTL, 

*/
#define MISSING NA_REAL

void outputResults(int nn, int iter, DATA* myData, QTL_INFO** all_qtls, double* resid, 
		   QTL_INFO** output_qtls, PRIORS* priors, MCMC_PARAM* myMCMC,
		   double* result,  int* output_idx)
{
  int iqtl,j;
  QTL_INFO* aqtl;
  double H, LOD;
  double residSS, val;
  static double logTwoPi = 1.837877066;
  int rj_mcmc = (myMCMC->revjump & MCMC_SELECT_PRIOR_FIELD);

  /*
   In following result was a double** to allow passage of new pointer to function mcmc.
 
  if (*result == NULL)
  {
     *allocRec = (int)(myMCMC->niter + myMCMC->burnIn) * ELEMENTS_PER_QTL_OUTPUT * MAXQTL/3;
     *result = Calloc(*allocRec, double);
  }
  if (*output_idx + MIN(1,myData->nQtl) > *allocRec)
  {
     Rprintf("Reallocating space\n");
     *allocRec = (int)(*allocRec * 1.25) ;
     *result = Realloc(result, *allocRec, double);  
  }
  */

  /* order QTLs for output */
  for(iqtl=1; iqtl<=myData->nQtl; iqtl++)
    {
      /* look for place to insert QTL in ordered sequence*/ 
      for(j=iqtl; ; j--)
	{
	  if (j==1 ||
	      all_qtls[iqtl]->chrom->num > output_qtls[j-1]->chrom->num ||
	      (all_qtls[iqtl]->chrom->num == output_qtls[j-1]->chrom->num &&
	       all_qtls[iqtl]->qptr->pos > output_qtls[j-1]->qptr->pos))
	    {
              output_qtls[j] = all_qtls[iqtl];
	      break;
	    }
	  else
	    output_qtls[j] = output_qtls[j-1];  /* open a space */
	}
    }

  if(myData->nQtl > 0) 
    {
      LOD = get_lod(nn, myData->sigmasq, myData->ybar, myData->y_var, resid);
      H = calc_h2(nn, myData->y, myData->mu, myData->sigmasq, resid);
    }
  else
    LOD = 0.0;

  if (!rj_mcmc)
    {
      myMCMC->idx++;
      residSS = calcResidSS(nn, resid);
      val = ((nn/2.0) * (logTwoPi + log(myData->sigmasq))) + 
	residSS/(2.0*myData->sigmasq);
      updateMean(myMCMC->idx, exp(val), &myMCMC->HM);

      val = priors->sig_a1 * log(1 + residSS/(2*priors->sig_a2));
      val += gammln2(priors->sig_a1 - (nn/2.0));
      val += ((nn/2.0)*(logTwoPi + log(priors->sig_a2)));
      val -= gammln2(priors->sig_a1);
      updateMean(myMCMC->idx, exp(val), &myMCMC->SHM);		
    }

  /* rewrite by Yandell in new format 30 mar 1999, PGA 11-18-2000 */
  if(myData->nQtl==0)
    {
/*      fprintf(write_res, "%d 0 1 . %15.4f %12.4f %12.4f ",
	      iter, LOD, myData->mu, myData->sigmasq);
      if (myMCMC->addParam == QTL_ADD_DOM || myMCMC->addParam == QTL_NONE)
	fprintf(write_res,"%12.4f %12.4f ", priors->var[QTL_ADD]/myData->y_var, 
		priors->var[QTL_DOM]/myData->y_var);
      else if (myMCMC->addParam == QTL_ADD)
	fprintf(write_res,"%12.4f . ", priors->var[QTL_ADD]/myData->y_var);
      else
	fprintf(write_res,". %12.4f ", priors->var[QTL_DOM]/myData->y_var);
      if (!rj_mcmc) fprintf(write_res,". . . 0 %12.4f %12.4f\n",myMCMC->HM, myMCMC->SHM);
      else fprintf(write_res,". . . 0\n");
*/
  /* write to result vector  PGA 8-21-2003 */
     result[(*output_idx)++] = iter;
     result[(*output_idx)++] = 0;
     result[(*output_idx)++] = 1;
     result[(*output_idx)++] = MISSING;
     result[(*output_idx)++] = LOD;
     result[(*output_idx)++] = myData->mu;
     result[(*output_idx)++] = myData->sigmasq;
     if (myMCMC->addParam == QTL_ADD_DOM || myMCMC->addParam == QTL_NONE)
       {
	 result[(*output_idx)++] = priors->var[QTL_ADD]/myData->y_var;
	 result[(*output_idx)++] = priors->var[QTL_DOM]/myData->y_var;
       }
     else if (myMCMC->addParam == QTL_ADD)
       {
	 result[(*output_idx)++] = priors->var[QTL_ADD]/myData->y_var;
	 result[(*output_idx)++] = MISSING;
       }
     else
       {
	 result[(*output_idx)++] = MISSING;
	 result[(*output_idx)++] = priors->var[QTL_DOM]/myData->y_var;
       }
     result[(*output_idx)++] = MISSING;
     result[(*output_idx)++] = MISSING;
     result[(*output_idx)++] = MISSING;
     result[(*output_idx)++] = 0;
    }
  else
    {
    for(iqtl=1; iqtl<=myData->nQtl; iqtl++)
      {
	aqtl = output_qtls[iqtl];
	
	/* output preamble ... including number of chromosome, LOD, mu, sigmasq */
	/*
	  if (iqtl == 1)
	  {
	  fprintf(write_res, "%d %d %d %d %15.4f %12.4f %12.4f", 
	  iter, myData->nQtl, iqtl, aqtl->chrom->num,
	  LOD, myData->mu, myData->sigmasq);
     	    if (myMCMC->addParam == QTL_ADD_DOM || myMCMC->addParam == QTL_NONE)
	      fprintf(write_res,"%12.4f %12.4f ", priors->var[ADD]/myData->y_var, 
		      priors->var[DOM]/myData->y_var);
	    else if (myMCMC->addParam == QTL_ADD)
	      fprintf(write_res,"%12.4f . ", priors->var[QTL_ADD]/myData->y_var);
	    else
	      fprintf(write_res,". %12.4f ", priors->var[QTL_DOM]/myData->y_var);
	  }
	else
	  fprintf(write_res, "%d %d %d %d . . . . . ", iter, myData->nQtl, 
		  iqtl, aqtl->chrom->num);

	if (aqtl->flag & QTL_ADD)
	  fprintf(write_res, " %12.4f",aqtl->a[QTL_ADD]);
	else 
	  fprintf(write_res, ". ");

	if (aqtl->flag & QTL_DOM)
	  fprintf(write_res, " %12.4f",aqtl->a[QTL_DOM]);
	else 
	  fprintf(write_res, " .");

	if (iqtl!=1)
	  {
	    if (!rj_mcmc) fprintf(write_res, " %10.4f . . .\n", aqtl->qptr->pos*100);
	    else fprintf(write_res, " %10.4f .  \n", aqtl->qptr->pos*100);
	  }
	else if (!rj_mcmc) 
	  fprintf(write_res, " %10.4f %8.4f %12.4f %12.4f\n",aqtl->qptr->pos*100,  H,
		  myMCMC->HM, myMCMC->SHM);
	else
	  fprintf(write_res, " %10.4f %8.4f\n",aqtl->qptr->pos*100,  H);
*/
	result[(*output_idx)++] = iter;
	result[(*output_idx)++] = myData->nQtl;
	result[(*output_idx)++] = iqtl;
	result[(*output_idx)++] = aqtl->chrom->num;
	if (iqtl==1)
	  {
	    result[(*output_idx)++] = LOD;
	    result[(*output_idx)++] = myData->mu;
	    result[(*output_idx)++] = myData->sigmasq;
	  }
	else
	  {
	    result[(*output_idx)++] = MISSING;
	    result[(*output_idx)++] = MISSING;
	    result[(*output_idx)++] = MISSING;
	  }
	if (myMCMC->addParam == QTL_ADD_DOM || myMCMC->addParam == QTL_NONE)
	  {
	    result[(*output_idx)++] = priors->var[QTL_ADD]/myData->y_var;
	    result[(*output_idx)++] = priors->var[QTL_DOM]/myData->y_var;
	  }
	else if (myMCMC->addParam == QTL_ADD)
	  {
	    result[(*output_idx)++] = priors->var[QTL_ADD]/myData->y_var;
	    result[(*output_idx)++] = MISSING;
	  }
	else
	 {
	   result[(*output_idx)++] = MISSING;
	   result[(*output_idx)++] = priors->var[QTL_DOM]/myData->y_var;
	 }
	
	if (aqtl->flag & QTL_ADD)
	  result[(*output_idx)++] = aqtl->a[QTL_ADD];
	else 
	  result[(*output_idx)++] = MISSING;
	
	if (aqtl->flag & QTL_DOM)
	  result[(*output_idx)++] = aqtl->a[QTL_DOM];
	else 
	  result[(*output_idx)++] = MISSING;
	
	result[(*output_idx)++] = aqtl->qptr->pos*100;
	if (iqtl==1)
	  result[(*output_idx)++] = H;
	else
	  result[(*output_idx)++] = MISSING;
      }
    }
}





double calc_h2(int nn, double* y, double mu, double sigmasq, double* resid)
{
  int i;
  double sum, mean, add;

  if (nn <=0) return 0;

  /* calculate the genetic variance, i.e. var(gene effects) */
  for (i=1, sum=0.0; i<= nn; i++) 
    sum += y[i] - mu - resid[i] ; 
  mean = sum/ nn;

  for (i=1, sum=0.0; i<= nn; i++) 
    {
      add = y[i] - mu - resid[i] - mean; 
      sum += add*add;
    }
  sum /= nn;       /* estimate for genetic variance */

  return sum / (sum + sigmasq);
	
}


void initVars(int nn, double** nmove, int* num_qtl, MCMC_PARAM* myMCMC,
	      double qtl_prior_mn, double cval)
{
  int i;

  for(i=1; i<=2; i++)
    nmove[C_BIRTH][i] = nmove[C_DEATH][i] = nmove[C_UPDATE][i] =
      nmove[C_EXCHANGE][i] = nmove[C_SWAP][i] = 0.0;

  for (i=0; i<=MAXQTL; i++)
    num_qtl[i] = 0;

}





void setupBurnIn (MCMC_PARAM* myMCMC, DATA* myData)
{
  double ratio, burnIn, preBurnIn;;

  
  /* BurnIn: a ratio of the number of iterations (10%) */
  if (myMCMC->burnIn <0) burnIn = BURN_IN * myMCMC->niter;
  else burnIn = myMCMC->burnIn * myMCMC->niter;

  if (myMCMC->preBurnIn <0) preBurnIn = myData->nChrom * 100;   
  else preBurnIn = myMCMC->preBurnIn * burnIn;

  if (myMCMC->burnIn < 0 || myMCMC->preBurnIn < 0)
    {
      /* PreBurnIn: get 50 birth proposals per chromosome on average */

      /* now check if PreBurnIn not long enough */
      ratio = preBurnIn/(burnIn * PRE_BURN_IN);
      if (ratio > 1) 
	{
	  burnIn *= ratio; 
	  myMCMC->niter = (int)(myMCMC->niter * ratio); 
	}
      else if (ratio < 1)
	{
	  preBurnIn = PRE_BURN_IN * burnIn;
	}
    }


  /*  printf("%lf %lf\n",burnIn, preBurnIn);
   */  
  myMCMC->burnIn = floor(burnIn);
  myMCMC->preBurnIn = floor(preBurnIn);
  /*  printf("%lf %lf\n",myMCMC->burnIn, myMCMC->preBurnIn);*/

}

/*
            h2 = LOD * log(10);
            h2 = 2*(h2-1)/(h2+myData->nn-1);
*/

void calcMeanVar(int nn, double* vals, double* mean, double* var)
{
  int i;
  double sum,t;

  for (i=1, sum=0.0; i<=nn; i++)
    sum += vals[i];
  *mean = sum/nn;

  for (i=1, sum=0.0; i<=nn; i++)
    {
      t = vals[i] - *mean;
      sum += t*t;
    }
  *var = sum/nn;
}


void updateMean(int i, double val, double* mean)
     /* i=1,2,.... 
   user must initialize mean=0, so when called with i=1, mean= val */
{
  double delta = (val - *mean) / i;
  *mean += delta;
}
