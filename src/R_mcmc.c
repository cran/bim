/************************************************************
 *
 * R_mcmc.c 
 * The C wrapper function for MCMC engine
 *
 * Written by: Hao Wu hao@jax.org, 08/15/2003
 #             Pat Gaffney   gamhna@ameritech.net  Added code from main.c 
 *
 ************************************************************/


#include "revjump.h"
#include "ranlib.h"

void R_mcmc(int *nind, int *nchr, /* number of individuals and chromsomes */
	    /* the following four items are for genetic map */
	    int *chrarray, int *nmararray, double *posarray, double *distarray,
	    int *genodata, /* genotype data */
	    int *nmar, /* number of markers per chromosome */
	    double *pheno, /* phenotype */
	    int *crosstype1, int *crosstype2, /* for cross type */
	    /* MCMC parameters */
	    double* burnin, double*preburnin, int*niter, int *nby,
	    int *seed, /* random number seed */
	    /* the following are priors */
	    double *meanadd, double *sdadd,
	    double *meandom, double *sddom,
	    double *meanmu, double *sdmu,
	    double *sigma2_a, double *sigma2_b,
	    char *priortype, double *priornqtl,
	    /* initial values for chain */
	    double *initmu, double *inits2,
	    /* return variable */
	    double *result)

{
  CHROMOSOME* chromInfo;
  DATA myData;
  WORK myWork;
  MCMC_PARAM myMCMC;
  PRIORS priors;
  params* theparams;
  int revjump;
  double* new_lambda=NULL;
  int* new_chrom=NULL;

   
  theparams = (params *)R_alloc(1, sizeof(params));
  theparams->cross = *crosstype1;
  theparams->crosst = *crosstype2;
  theparams->chrom = *nchr;
  theparams->error = NULL;

 
  if (!strcmp(priortype,"geometric"))
     revjump=SELECT_NQTL_GEOMETRIC;
  else if (!strcmp(priortype,"uniform"))
     revjump=SELECT_NQTL_UNIFORM;  
  else if (!strcmp(priortype,"fisch"))
     revjump=SELECT_NQTL_FISCH;
  else if (!strcmp(priortype,"poisson"))
     revjump=SELECT_NQTL_POISSON;
  else 
  {
     revjump=SELECT_NQTL_GEOMETRIC;
  }


  /* setup number generator */
  while (*seed < 0) *seed += 2147483563;
  if (*seed == 0) *seed = time(NULL);    
  if ((*seed&16) == 0) setall(*seed/4+13, *seed/3+211);
  else setall(*seed/3+211, *seed/4+13);
  advnst((*seed & 31)+3);   
  /*  theparams->seed = *seed; */


  setupMCMC(&myMCMC, &myData, &myWork,  &priors, 
	    theparams, genodata, 
	    chrarray, nmararray, posarray, distarray,
	    nmar, pheno, &chromInfo, &revjump, nind,
	    burnin, preburnin, niter, nby,
	    initmu,inits2,meanmu,sdmu,sigma2_a,sigma2_b,
	    meanadd, sdadd, meandom, sddom, priornqtl);
  
  
  mcmc (&myMCMC, &myData, &priors, chromInfo, &myWork, result);
}


void setupMCMC(MCMC_PARAM* myMCMC, DATA* myData, WORK* myWork, PRIORS* priors,
	       params* theparams, int* genodata,
	       int *chrarray, int *nmararray, double *posarray, double *distarray,
	       int* markersPerChrom, double* pheno, CHROMOSOME** pChromInfo,
	       int* revjump,int* nind, double* burn_in,
	       double* preburn_in,int* niter,int* nby, double* mu, double* sigmasq,
	       double* meanmu, double* sdmu, double* siga1, double* siga2,
	       double* meanadd, double *sdadd,double* meandom, double*sddom,
	       double* qtlmean)
{
  double* new_lambda=NULL;
  int* new_chrom=NULL;
  
  /* a little help optimizing code for 2-genotype crosses, like BC, DH */
    switch (theparams->cross) {
    case 1:  myData->gmiss = -1; break;
    case 2:  myData->gmiss =  1; break;
    case 5:  myData->gmiss =  0;  break;
    default: myData->gmiss = -2;
    };
    
   setupTraitData(myData, *nind, pheno);
   setupChromosomes(myData, myMCMC, pChromInfo, theparams, genodata, 
	                chrarray, nmararray, posarray, distarray, markersPerChrom);
   setupWork(myData->nn, myWork);

   /* miscellaneous */
   myData->theparams = theparams;
   myData->bp=myData->dp=myData->prior_ratio=(double*)NULL;  /* initialize */
   
  
   getdata(myMCMC, myData, priors, *pChromInfo, myWork,revjump,
	       burn_in,preburn_in,niter,nby,mu,sigmasq,meanmu,sdmu,siga1,siga2,
		   meanadd, sdadd, meandom, sddom, qtlmean ); 

   setupQtl(myData, myMCMC,  myWork, *pChromInfo, priors, new_chrom, new_lambda,
 		    theparams);
   

   priors->priorDistribution= myMCMC->revjump & 0xF;
   
   /* setup ... */
   if ((myMCMC->revjump & MCMC_SELECT_MIX) == MCMC_OPTIMIZE_80_20) 
     myMCMC->cval = -0.8;
   else if ((myMCMC->revjump & MCMC_SELECT_MIX) == MCMC_OPTIMIZE_70_30) 
     myMCMC->cval = -0.7;
   else if ((myMCMC->revjump & MCMC_SELECT_MIX) == MCMC_OPTIMIZE_50_50) 
     myMCMC->cval = -0.5;
   else 
     myMCMC->cval = -0.9;  /* default */

   /* decide whether to make the birth/death (BD) probs, bp+dp, exactly
      equal to myMCMC->cval (default) or to make it an inequality */
   if (myMCMC->revjump & MCMC_SELECT_BD_NOTEXACT) myMCMC->cval = -myMCMC->cval;
   
   setupBirthDeathProbs(myMCMC->revjump, priors->qtl_mean, 
			myMCMC->cval, &myData->bp, &myData->dp, &myData->prior_ratio);/* here*/


   myMCMC->HM = myMCMC->SHM = 0;
   myMCMC->idx = 0;
   

   if (myMCMC->revjump & MCMC_SELECT_SIMPLE_GIBBS) 
     {
       myMCMC->revjump |= MCMC_SELECT_MOVE_NONE;
       myMCMC->revjump |= SELECT_SAMPLE_FISCH;
   }
}




void setupWork(int nn, WORK* myWork)
{
  int i;
  
  
   /* residuals */
  myWork->resid = dvector(1,nn);
  myWork->newResid = dvector(1,nn);
  
  
  myWork->u = dvector(1,MAX_CHOL);
  myWork->pmean = dvector(1,MAX_CHOL);
  myWork->pvar = dvector(1,MAX_CHOL);
  myWork->p = dvector(1,MAX_CHOL);
  myWork->XtY = dvector(1,MAX_CHOL);
  myWork->new_XtY = dvector(1,MAX_CHOL);
  myWork->XtX = dmatrix(1,MAX_CHOL, 1, MAX_CHOL);
  myWork->new_XtX = dmatrix(1,MAX_CHOL, 1, MAX_CHOL);
  myWork->chol = dmatrix(1,MAX_CHOL, 1, MAX_CHOL);
  myWork->perm_num = ivector(1,MAXQTL);
  
  /* diagnostics */
  myWork->minProbs = dvector(1,nn);
  myWork->avgMinProbs = dvector(1,nn);
  
  /* for qtl updating */
  myWork->new_r = dvector(-1,1);
  myWork->log_norm_const = dvector(1,nn);
  myWork->output_qtls = (QTL_INFO**)R_alloc(MAXQTL, sizeof(QTL_INFO*));
  myWork->output_qtls--;
  for (i=1; i<= MAXQTL; i++) myWork->output_qtls[i]=(QTL_INFO*)NULL;
  
  myWork->mod_effect = dvector(0,MAX_CHOL);
  myWork->weight = dvector(0,MAX_CHOL);
  myWork->work = dvector(0,MAX_CHOL);
  myWork->oldWeight = dvector(0,MAX_CHOL);
  myWork->oldVar = dvector(0,MAX_CHOL);
  myWork->weight[0] = 1.0;          /* the weight for the prior mean variance */
}




void noramlizeBirthDeath(double* bp, double* dp, double cval)
  /* if cval is positive, then we ensure that they are no more
   * than cval.   The latter ensures that the prior ratio cancels
   * out the ratio of death to birth proposal probabilities in the
   * acceptance ratio for RJ-MCMC (and is standard Green 1995).
   *
   * If cval is negative, we ensure that death and birth probs 
   * add up to exactly cval 
   */
{
    int i;
    double mval;
    
    if (cval > 0 && cval <= 1)
      {
	mval = 0;
	for (i=0; i<=MAXQTL; i++)
	  mval = MAX(mval, (bp[i] + dp[i]));
	mval = cval/mval;
	
	for (i=0; i<=MAXQTL; i++)
	  {
	    bp[i] *= mval; 
	    dp[i] *= mval;
	  }
      }
    else {
      cval = -cval;
      if (cval <= 0) cval = 0.9;
      for (i=0; i<=MAXQTL; i++)
	{
	  bp[i] = cval * bp[i]/(bp[i] + dp[i]);
	  dp[i] = cval - bp[i];
	}
    }
}




void setupBirthDeathProbs(int revjump, double qtl_param, 
			  double cval,
			  double** bp, double** dp, double** priorRatio)
{
  int i;
  int priorType = (revjump & MCMC_SELECT_PRIOR_FIELD);
  int proposalType = (revjump & MCMC_SELECT_PROPOSE_FIELD);
  
  /* birth/death proposal probabilities */
   if (!*bp) *bp = dvector(0,MAXQTL);
   if (!*dp) *dp = dvector(0,MAXQTL);
   if (!*priorRatio) *priorRatio = dvector(1,MAXQTL);
   (*dp)[0] = (*bp)[MAXQTL] = 0.0;
   
   if (priorType==0)
     {
       for (i=0; i<=MAXQTL; i++) (*bp)[i] = (*dp)[i] = (*priorRatio)[i] = 0.0;
       return;
     }
   
   for (i=0; i<MAXQTL; i++)
     {
       switch (priorType) 
	 {
	 case SELECT_NQTL_POISSON:
		  (*priorRatio)[i+1] = qtl_param/(i+1);
		  break;
		  
	 case SELECT_NQTL_GEOMETRIC:
	   (*priorRatio)[i+1] = 1.0/qtl_param;
  	      break;
	      
	 case SELECT_NQTL_FISCH:
	   (*priorRatio)[i+1] = qtl_param/(i+1)/(i+1);
	   break;
	   
	 case SELECT_NQTL_UNIFORM:
	   (*priorRatio)[i+1] = 1;
	   break;
	 default: Rprintf("Should have poisson,geometric,"
			 "fisch or uniform for prior type. "
			 "Choosing poisson as prior.\n");
	   (*priorRatio)[i+1] = qtl_param/(i+1);
	 };
       
       switch (proposalType)		  
	 {
	 case SELECT_PROPOSE_PRIOR:
	   (*bp)[i] = MIN(1.0, (*priorRatio)[i+1]);
	   (*dp)[i+1] = MIN(1.0, 1.0/(*priorRatio)[i+1]);
	   break;
	 case SELECT_PROPOSE_POISSON:
	   (*bp)[i] = MIN(1.0, qtl_param/(i+1));
	   (*dp)[i+1] = MIN(1.0, (i+1)/qtl_param);
	   break;
	 case SELECT_PROPOSE_GEOMETRIC:
	   (*bp)[i] = MIN(1.0, 1.0/qtl_param);
	   (*dp)[i+1] = MIN(1.0, qtl_param);
	   break;
	 case SELECT_PROPOSE_FISCH:
	   (*bp)[i] = MIN(1.0, qtl_param/(i+1)/(i+1));
	   (*dp)[i+1] = MIN(1.0, (i+1)*(i+1)/qtl_param);
	   break;
	 case SELECT_PROPOSE_UNIFORM:
	   (*bp)[i] = 1;
	   (*dp)[i+1] = 1;
	   break;
	 default: Rprintf("Should have poisson,geometric,"
			 "fisch or uniform. "
			 "Choosing poisson default.\n");
	   (*bp)[i] = MIN(1.0, (*priorRatio)[i+1]);
	   (*dp)[i+1] = MIN(1.0, 1.0/(*priorRatio)[i+1]);
	  };
       
       
       
     }
   noramlizeBirthDeath(*bp, *dp, cval);
}










void setupQtl(DATA* myData, MCMC_PARAM* myMCMC, WORK* myWork, 
	      CHROMOSOME* chromInfo, PRIORS* priors, 
	      int* new_chrom, double* new_lambda, 
	      params* theparams)			  
{
  QTL_INFO* newQtl;
  CHROMOSOME* chrom;
  int nn = myData->nn;
  int i, c, lmark;
  
  /* qtl info .. we allocate two more than normal for work variables */
  myData->myQtls = (QTL_INFO**)R_alloc(MAXQTL+2, sizeof(QTL_INFO*));
  for (i=0; i< MAXQTL+2; i++) myData->myQtls[i]=(QTL_INFO*)NULL;
  createQtl(myData->nn, 0, &myData->myQtls[0], NULL, 0, 0.0, QTL_NONE, NULL, NULL);
  
  
  myData->na[ADD]=myData->na[DOM]=0;
  
  for( i=1; i<=myData->nQtl; i++ ) 
    {
      /* allow for errors, and special case 
	 new_chrom[i] == 0 => use qtl index i as chromosome number
	 new_lambda[i] == 0 => use middle of chromosome as starting position
      */
      c = new_chrom[i];
      if (c < 1 || c > myData->nChrom)
  	  c = ignuin(1,myData->nChrom);
      chrom = &chromInfo[c];
      
      if (new_lambda[i] <= 0.0 || new_lambda[i] >= chrom->chromLen)
	get_local_locus(chrom, &lmark, &new_lambda[i]);
      else
	lmark = binSearch(chrom->nMark, chrom->mark_pos, new_lambda[i]);
      
      newQtl = createQtl(nn, i, &myData->myQtls[i], chrom, 
			 lmark, new_lambda[i], myMCMC->addParam, NULL, NULL);
      if (newQtl->flag & QTL_ADD) myData->na[QTL_ADD]++;
      if (newQtl->flag & QTL_DOM) myData->na[QTL_DOM]++;
      
      
      if (myMCMC->addParam == QTL_NONE) {newQtl->flag = QTL_ADD_DOM;}  
      /* force it to fully fit*/
      
      initQtl(myData->nn,  newQtl, 
	      theparams, myData->bc, myData->gmiss, myMCMC->revjump);
      newQtl->chrom->nQtl++;
    }
  
  if (myData->nQtl > 0) free_dvector(new_lambda,1,myData->nQtl);
  
  update_effects(nn, myData->nQtl, myMCMC->revjump, myData->y, 
		 &myData->mu, myData->sigmasq, myData->na, 
		 myData->myQtls, priors, myWork);
  
}





void setupChromosomes(DATA* myData, MCMC_PARAM* myMCMC, CHROMOSOME** pChromInfo,
		      params* theparams, int* genodata, 
		      int *chrarray, int *nmararray, double *posarray, 
		      double *distarray, int* markersPerChrom)
{
  int i,c, idx, mark_idx, k;
  mygenome *mystartptr, *prev;
  CHROMOSOME* chromInfo;
  CHROMOSOME* chrom;
  int nn = myData->nn;
  
  /* setup the number of chromosomes to analyze */
  myData->nChrom = theparams->chrom;
  myMCMC->offset = 1;   
					        
  /* allocate space for chromosome information */
  *pChromInfo = (CHROMOSOME*)R_alloc(myData->nChrom, sizeof(CHROMOSOME)); 
  chromInfo = --(*pChromInfo);   /* so indexing is 1,2,.... */
  myData->myChroms = chromInfo;
  myData->totalMark = 0;
  myData->totalChromLen = 0;
  myData->chrom_pos = dvector(1,myData->nChrom+1);
  myData->chrom_pos[1] = 0.0;
  
  
  
  idx=0;
  mark_idx=0;
  markersPerChrom--;
  for (c=1; c<= myData->nChrom; c++) 
    {
      chrom = &chromInfo[c];
      chrom->nMark = markersPerChrom[c];
      chrom->nQtl=0;
      chrom->num=-1;
	  
      chrom->mark_genome = (mygenome**)R_alloc(chrom->nMark+1, sizeof(mygenome*));
      chrom->mark_pos = dvector(1,chrom->nMark);
      chrom->mark_genome--;
      prev = (mygenome*)NULL;
      
      /* make a reference array for marker info */
      for (k=1; k<= chrom->nMark; k++) 
	{
	  mystartptr = chrom->mark_genome[k] = (mygenome*)R_alloc(1,sizeof(mygenome));
	  if (chrom->num <0)
	    chrom->num = chrarray[idx];
	  else if (chrom->num != chrarray[idx])
	    Rprintf("Error -- reading map info for chromosome %d and found entry for chromosome %d \n",
		   chrom->num, chrarray[idx]);
	  mystartptr->chrom = chrom->num;
	  mystartptr->markr = nmararray[idx];
	  mystartptr->pos = posarray[idx]/100;
	  mystartptr->dist = distarray[idx]/100;
	  idx++;
	  mystartptr->next=NULL;
	  
	  
	  mystartptr->prev = prev;
	  if (prev) prev->next = mystartptr;
	  
	  
	  chrom->mark_pos[k] = mystartptr->pos;
	  prev = mystartptr;
	  
	  chrom->mark_genome[k]->genotype = &genodata[nn*mark_idx];
	  chrom->mark_genome[k]->genotype--;
	  mark_idx++;
	}
      
      if (mystartptr->dist>0) 
	{
	  mystartptr = chrom->mark_genome[k] = (mygenome*)R_alloc(1, sizeof(mygenome));
	  mystartptr->chrom = prev->chrom;
	  mystartptr->markr=0;
	  mystartptr->dist=0.0;
	  mystartptr->pos= prev->pos + prev->dist;
	  mystartptr->prev=prev;
	  prev->next = mystartptr;         
	}
      
      
      /* capture some overall statistics (why not!) */
      /* keep some chromosome length statistics (needed for position proposal) */
      chrom->chromLen = mystartptr->pos;
      myData->chrom_pos[c+1] = myData->chrom_pos[c] + chrom->chromLen;
      
      myData->totalMark += chrom->nMark;
      myData->totalChromLen += chrom->chromLen;
      
      /* set up record of chromosomes on the QTL */
      chrom->qtls = (QTL_INFO**)R_alloc(MAXQTL, sizeof(QTL_INFO*));	 
      chrom->qtls--;
      for (i=1; i<= MAXQTL; i++) chrom->qtls[i]=(QTL_INFO*)NULL;
    }
}




void setupTraitData(DATA* myData, int nn, double* pheno)
{  
  myData->nn = nn;
  myData->y = pheno-1;
  calcMeanVar(myData->nn, myData->y, &myData->ybar, &myData->y_var);
  myData->y_stdev = sqrt(myData->y_var);
}  
