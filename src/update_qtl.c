/*****************************************************
  File:         update_qtl.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:
    Update qtl genotypes from Bernoulli probability

*****************************************************/

/*-------------------------------------------------------------------------
 *
 *
 *
 *
 *-------------------------------------------------------------------------*/


#include "revjump.h"
#include "ranlib.h"

void bselfed_f_tpm(double (*sftpm)[3], double r, int t, int bc);
void selfed_f_tpm(double (*sftpm)[3], double r, int t);
void assign_q(double *ql,params *theparams);
void a_Mv(double *a, double (*M)[3], double *v);

void AdotB2(double **rr, double (*aa)[3], double (*bb)[3]);
void AinR2(double (*rr)[3], double **aa);
void bselfed_f_tpm2(double **sftpm, double r, int t, int bc);
void selfed_f_tpm2(double** sftpm, double r, int t);


double update_lambda_qtl(QTL_INFO* aqtl, QTL_INFO** all_qtls, DATA* myData, 
						 WORK* myWork, int revjump)
{
    int newMark;
	double newPos;
    double accept_sum=0;
	double accept_lam;


		setValidFlag(aqtl, revjump);        /* invalidate the cache of log probabilities (of
								      QTL given flanking marker) for all individuals
									  if warranted */

        get_new_pos(aqtl->chrom, aqtl->qptr->pos, 0.03, &newMark, &newPos);	     

        accept_lam = accept_new_lambda(aqtl, all_qtls, myData->nn, 
			                           newMark, newPos, myData->gmiss,
		   					           myData->theparams, myData->bc, revjump);


#ifndef MCMC_OMIT_DEBUG
		checkIntegrity(myData->nQtl, aqtl->chrom);
#endif

        Gibbs_update_geno(aqtl, myData->nn, &myWork->resid,  myData->sigmasq,
			              myData->theparams, myData->bc, myData->gmiss, 
						  myWork->new_r, p_igenotype(aqtl), 
						  &myWork->newResid, myWork->log_norm_const, revjump);			 


		accept_sum += accept_lam/myData->nQtl;


#ifndef MCMC_OMIT_DEBUG
	    checkResid(myData->nn, myData->nQtl, myData->mu, myData->y, myData->myQtls, myWork->resid,
	               myData->gmiss);
#endif


	
	
	return (double)accept_sum;
}



void setValidFlag(QTL_INFO* aqtl, int revjump)
/* We have two lookup structures:

     (1) a matrix (3x3x3) giving the probability for any QTL genotype given flanking markers, 
	     condProb
	 (2) a vector (nn) giving probability for current QTL genotype for each individual given
	     all marker information (on chromosome), prob.

     (1) remains correct unless the flanking marker is a QTL which has just moved.  The vector
	 of probabilities in (2) is not as simple, since it's main use is for the mising data case
	 (to save recomputing probabilities using QTL Cartographer's calc_cond_prob routine).
	 One sufficent condition for validity is that there is only one QTL on the chromosome (a
	 frequent occurrence, so worth coding for).

     (2) is checked in setupTable
*/
{
	aqtl->probValid = aqtl->probValid && (aqtl->chrom->nQtl==1) &&
		              !(revjump & MCMC_OPTIMIZE_NONE);
}


void get_new_pos(CHROMOSOME* chrom, double oldPos, double delta, int* newMark, double* newPos)
/* all distances in Morgans */
{
  if (delta > chrom->chromLen) delta = chrom->chromLen/2.0;

  *newPos = oldPos + genunf(-1.0,1.0) * delta;
  if (*newPos < 0) *newPos = -*newPos;
  if (*newPos > chrom->chromLen) *newPos = 2 * chrom->chromLen - (*newPos);
  *newMark = binSearch(chrom->nMark, chrom->mark_pos, *newPos);
}


double get_new_lambda(int nMark, double lambda)
{
  double low, high, propose;

    low = lambda - 1.0;
    high = lambda + 1.0;

	/* if the value proposed lies outside the bounds, we "reflect" the value,
	   i.e. if propose is 0.8, we set it equal to 1.2 (equal distance from our
	   minimum value of 1.0.  
	   In MCMC, this ensures that the ratio of proposal densities is 1.  Note 
	   that it does affect the density for proposal, but this is not important 
	   (only the ratio).
	*/
	propose = low + (high-low)*genunf(0, 1);

	while (propose < 1.0 || propose > nMark)
	{
	  if( propose < 1.0 ) 
	     propose = 2-propose;
      else if( propose > nMark ) 
  	    propose = nMark + nMark-propose;
	}

	return propose;
}









int setNewFlank(mygenome* qptr, int* prevRenew, int* nextRenew,
				double* prevDist, double* nextDist, int revjump)
/*-------------------------------------------------------------------------
 *
 *  determines the flanking marker (or qtl) genotypes to left and right
 *  of specified qtl (aqtl).  Stores them in aqtl->flankLeft and 
 *  aqtl->flankRight respectively.  
 *
 *  Finally, determines if these flanking markers have changed, and if so, 
 *  the recombination rates between flanking markers and QTL (aqtl->rLeft
 *  and aqtl->rRight) and between the markers (aqtl->R) is computed.
 *
 *  If you do not want recombination rate to be calculated (between
 *  QTL and flanking marker/QTLs, set aqtl->flankRight and aqtl->flankLeft
 *  to (int*)NULL.  When proposing a new QTL do not do this.  This is present
 *  simply to check if a flanking marker was a QTL and it has changed 
 *  position.
 * 
 *-------------------------------------------------------------------------*/
{
   mygenome* prev, *next;


   prev = qptr->prev;  /* note we guarantee a prev/next marker */
   next = qptr->next;  /* it may have missing genotype data    */

   *prevRenew = (*prevDist != prev->dist);
   *nextRenew = (*nextDist != qptr->dist);
   *prevDist =  prev->dist;
   *nextDist = qptr->dist;
/*   
   safer code 
   ==========
   *prevRenew = (*prevDist < 0) || (prev->markr < 0);
   *nextRenew = (*nextDist < 0) || (next->markr < 0);
*/
   if ((revjump & MCMC_SELECT_OPTIMIZE) == MCMC_OPTIMIZE_NONE) 
       return 0;     /* since we won't ever use the table, prevent generate */	                   
   else if (revjump & MCMC_OPTIMIZE_NOTABLE) 
       return 1;  /* force recalculation it from being computed */

   if (!(*prevRenew) && !(*nextRenew)) return 0;


   return 1;
}


double calc_flank_dist(mygenome* qptr)
{
	return qptr->dist + qptr->prev->dist;
}



  





void setLambda(QTL_INFO* aqtl, int mark, double pos, int qtlNum)
/*-------------------------------------------------------------------------
 *
 *  
 *  qtlNum   > 0 represents index of ordered QTL along chromosome
 *  pos      position of QTL
 *
 *
 *-------------------------------------------------------------------------*/
{
    if (aqtl->chrom) insertQtl(aqtl->qptr, mark, pos, aqtl->chrom, qtlNum);
}



void setUpGeno(mygenome* qptr, int** prevGeno, int** nextGeno, int** geno)
{
	*prevGeno = &(qptr->prev->genotype[1]);
	*nextGeno = &(qptr->next->genotype[1]);
	*geno = &(qptr->genotype[1]);
}



void SwapGenome(mygenome** g1, mygenome** g2)
{
	mygenome* g = *g1;
	*g1=*g2;
	*g2=g;
}


void Swap3DTable(double**** g1, double**** g2)
{
	double*** g = *g1;
	*g1=*g2;
	*g2=g;
}

void Swap2DTable(double*** g1, double*** g2)
{
	double** g = *g1;
	*g1=*g2;
	*g2=g;
}

void SwapIVec(int** g1, int** g2)
{
	int* g = *g1;
	*g1=*g2;
	*g2=g;
}


void SwapDVec(double** g1, double** g2)
{
	double* g = *g1;
	*g1=*g2;
	*g2=g;
}


void SwapInt(int* g1, int* g2)
{
	int g = *g1;
	*g1=*g2;
	*g2=g;
}

void SwapDble(double* g1, double* g2)
{
	double g = *g1;
	*g1=*g2;
	*g2=g;
}

void swapQtl(QTL_INFO** g1, QTL_INFO** g2)
{
	QTL_INFO* g = *g1;
	*g1=*g2;
	*g2=g;
	SwapInt(&(*g1)->qptr->markr, &(*g2)->qptr->markr);
}




void copySwapGenoTo(QTL_INFO* newQtl, QTL_INFO* aqtl)
/*-------------------------------------------------------------------------
 *
 *  this copies minimal information.  Must call setLambda to insert
 *  the qtls's qptr record into the genome list.
 *
 *
 *  WARNING:  qptr genotype records are swapped.  This must be reset
 *            if the jump is not taken.
 *
 *
 *-------------------------------------------------------------------------*/
{
  newQtl->chrom = aqtl->chrom;
  newQtl->qptr->chrom = aqtl->qptr->chrom;
  newQtl->qptr->markr = aqtl->qptr->markr;
  newQtl->a[QTL_ADD] = aqtl->a[QTL_ADD];
  newQtl->a[QTL_DOM] = aqtl->a[QTL_DOM];
  newQtl->w[QTL_ADD] = aqtl->w[QTL_ADD];
  newQtl->w[QTL_DOM] = aqtl->w[QTL_DOM];
  newQtl->flag = aqtl->flag;

  SwapIVec(&newQtl->qptr->genotype, &aqtl->qptr->genotype);

  /* in ate entries */
  newQtl->nextDist = -1;
  newQtl->prevDist = -1;
  newQtl->qptr->pos = -1;
  newQtl->probValid = 0;
}


void swapQtlData(QTL_INFO* Qtl1, QTL_INFO* Qtl2)
{
	/* swap the chromosome record */
    CHROMOSOME* c = Qtl1->chrom;
	Qtl1->chrom = Qtl2->chrom;
	Qtl2->chrom = c;

    SwapGenome(&Qtl1->qptr, &Qtl2->qptr);    /* note that genome record contains genotype */
	Swap3DTable(&Qtl1->log_condProb, &Qtl2->log_condProb);    
	Swap2DTable(&Qtl1->transProbL, &Qtl2->transProbL);    
	Swap2DTable(&Qtl1->transProbR, &Qtl2->transProbR);    
	SwapDVec(&Qtl1->a, &Qtl2->a);
	SwapDVec(&Qtl1->w, &Qtl2->w);
	SwapInt(&Qtl1->flag, &Qtl2->flag);
	SwapInt(&Qtl1->nParam, &Qtl2->nParam);
	SwapDble(&Qtl1->prevDist, &Qtl2->prevDist);    
	SwapDble(&Qtl1->nextDist, &Qtl2->nextDist);    
	Swap2DTable(&Qtl1->log_prob, &Qtl2->log_prob);    
	Swap2DTable(&Qtl1->missing_prob, &Qtl2->missing_prob);    
	SwapInt(&Qtl1->probValid, &Qtl2->probValid);
}






/*--------------------------------------------------------------------------
 *  Function accept_new_lambda
 *
 *  essentially get the acceptance for the ith qtl
 *  on the given chromosome, moving from its current
 *  position to new_lambda[iqtl].
 *
 *  gmiss is -1(aa), 0(Aa) and 1(AA).
 *
 *--------------------------------------------------------------------------
 */

int accept_new_lambda(QTL_INFO* aqtl, QTL_INFO** all_qtls, int nn, 
					  int newMark, double newPos, int gmiss,
                      params* theparams, int bc, int revjump)
{
  double log_accept_prob;
  double uni_ran;
  int i;
  double** newProb, **oldProb;
  int* geno = &aqtl->qptr->genotype[1];
  QTL_INFO* newQtl;


#ifdef CHECK_TABLE
  double pi[3];
  int ii;
#endif

  mygenome* qptr = aqtl->qptr;
  int qtlNum = -qptr->markr;        /* number of marker being updated */

  /* calculate prob of QTL genos given flanking markers */
  genProbs(nn, theparams, bc, gmiss, aqtl, revjump);

  /* now temporarily replace the existing one */
  removeQtl(qptr);  

  /* setup new Qtl */
  newQtl = createQtl(nn, -qptr->markr, &all_qtls[0], aqtl->chrom, 
	                 newMark, newPos, aqtl->flag, aqtl->a, aqtl->w);	                 

  /* a clever pointer swap, but must restore later if move accepted (via
     swapQtlData) or not (via SwapIVec) */
  SwapIVec(&(newQtl->qptr->genotype), &(aqtl->qptr->genotype));  

  /* calculate prob of QTL genos given new flanking markers */
  genProbs(nn, theparams, bc, gmiss, newQtl, revjump);


  /* ratio of priors */
  log_accept_prob = 0.0;

  oldProb = aqtl->log_prob;
  newProb = newQtl->log_prob;

#ifdef CHECK_TABLE
  ii = ignuin(1,nn);
  calc_cond_prob2(theparams, bc, newQtl->qptr, 
	              ii, &pi[2], &pi[1], &pi[0]);
  for (i=0; i<=2; i++) 
	  if (i-1 != gmiss && pi[i] > 1e-10 && 
	      !EQUALS((log(pi[i])-newProb[ii][i-1])/newProb[ii][i-1],0.0))
  {
	  printf("Difference in log_prob in accept_new_lambda %20.17lf\n", 
		      log(pi[i])- newProb[ii][i-1]);
	  printf("Optimization failing, rerun program with optimization disabled\n");
	  printf("Set first line in nval.dat to 0x30000001 if it was 0x00000001,\n");
	  printf("or 0x30000002 if it was 0x00000002,etc.\n");
	  exit(1);
  }
#endif

  
  /* finally calculate the ratio */
   for(i=1;i<=nn;i++, geno++)
      log_accept_prob += (newProb[i][*geno] - oldProb[i][*geno]);
 
  uni_ran = genunf(0.0, 1.0);
  uni_ran = log(uni_ran);

  if( uni_ran < log_accept_prob ) 
  {
	  swapQtlData(aqtl, newQtl);
	  return 1;
  }
  else
  {
     /* a clever pointer swap ... saves a copy of nn items earlier */
     SwapIVec(&(newQtl->qptr->genotype), &(aqtl->qptr->genotype));
  }
  
  
  /* restore it to the way it was  remove inserted qtl genome record */
  /* put back old record */
  removeQtl(newQtl->qptr);           
  restoreQtl(qptr);            
  return 0;
}







void Gibbs_update_geno(QTL_INFO* aqtl, int nn, 
			 	  	   double** p_resid, double sigmasq,  
                       params* theparams, int bc, int gmiss,
					   double* new_r, 
					   int** p_genotype, double** p_newResid,
					   double* log_norm_const, int revjump)
{
/* propose the qtl genotype, then accept new genotype and residuals
   (i.e. swap existing ones with new ones)
*/
	propose_qtl_geno(aqtl, nn, 
			 		 *p_resid, sigmasq,  
                     theparams, bc, gmiss,
					 new_r, *p_genotype, *p_newResid, log_norm_const, revjump);

	SwapDVec(p_resid, p_newResid);
	SwapIVec(p_igenotype(aqtl), p_genotype);
}




void propose_qtl_geno(QTL_INFO* aqtl, int nn, 
			 		 double* resid, double sigmasq,  
                     params* theparams, int bc, int gmiss,
					 double* new_r, int* genotype, double* newResid,
					 double* log_norm_const, int revjump)
/* 
   proposes new values for the genotype for aqtl.  Updates the values in
   genotype, and provides new residuals in newResid.  
   
	 The caller must decide whether to accept or reject these proposals,
   and update (using SwapDVec) the old with new residuals, and the old 
   (aqtl->qptr->genotype) with the new genotypes.  In the case Gibbs 
   sampler, this occurs with 100% probability.
*/
{
  int i;
  int* geno, *newgeno;
  double log_qprob[3]; 

  genProbs(nn, theparams, bc, gmiss, aqtl, revjump);

  geno = &(aqtl->qptr->genotype[1]);
  newgeno = &genotype[1]; 


  for(i=1; i<=nn; i++, geno++, newgeno++)
  {
    log_norm_const[i] = gen_qprob(aqtl->a, *geno, sigmasq, resid[i], 
	                               aqtl->log_prob[i], &log_qprob[1], new_r, gmiss);

    GenGenotype(gmiss, &log_qprob[1], newgeno);
	newResid[i] = new_r[*newgeno]; 			
  }
}






/******************************************
        Function gen_qprob
******************************************/

double gen_qprob (double* a, int geno, double sigmasq,             
			  double resid,
              double* log_prob, double *log_qprob, double* new_r, int gmiss)
{
/*
   computes the conditional probability of Qtl genotype given marker and trait data
   and returns the log of the normailizing constant
*/
  int igen;
  double denominator=0.0;
  double log_norm_const;


  for(igen=-1;igen<=1;igen++) 
	  if (igen != gmiss) new_r[igen] = resid + EFFECT_VALUE(geno, a) - EFFECT_VALUE(igen, a);

  for(igen=-1; igen<=1; igen++) if (igen != gmiss)
  {
     log_qprob[igen] = log_prob[igen]  - 0.5*new_r[igen]*new_r[igen]/sigmasq;
     denominator += exp(log_qprob[igen]);
  }

  log_norm_const = -log(denominator);

  /* normalize, so that sum of these probabilities is one */
  for(igen=-1;igen<=1;igen++)
	  if (igen != gmiss) log_qprob[igen] += log_norm_const;
        


  return log_norm_const;
}











void setupTable(int gmiss, params* theparams, 
				int renewPrev, int renewNext, QTL_INFO* aqtl)
/*-------------------------------------------------------------------------
 *
 *  Generates conditional probability log_condProb[i][j][k] for QTL with 
 *  genotype k, given markers lRec to left with genotype i and rRec to 
 *  right with genotype j (i,j,k are -1(aa) to 1(AA).
 *
 *
 *  gmiss can be -1(aa), 0(Aa) and 1(AA).
 *  setNewFlank MUST be called
 *
 *  note that MissMark.c makes 0=>AA, 1=>Aa and 2=>aa   
 *  also transProbL[i][j] = Prob(X=j | Y=i) where X and 
 *  Y and adjacent loci, hence [1-i] when one would     
 *  naturally expect [i] (to convert 2 to -1,.., 0 to 1 
 *-------------------------------------------------------------------------*/
{
	int i,j,k;
	double p, denom;
	double* pr;
	double** transProbL = aqtl->transProbL;
	double** transProbR = aqtl->transProbR;
	double*** log_condProb = aqtl->log_condProb;

	double lRec = haldane(aqtl->prevDist);
	double rRec = haldane(aqtl->nextDist);
 
   
   if (renewPrev) calc_trans_prob(theparams, lRec, transProbL);
   if (renewNext) calc_trans_prob(theparams, rRec, transProbR);

	gmiss = 1-gmiss;  /* coding AA=0 (in calc_trans_prob which is
	                     almost a direct copy of QTLCart's markov_rec)
	                     versus AA=1 in our (and remainder of QTLCart's)
	                     code.  We convert from 1 to 0, 0 to 1, and -1 to 2
	                  */


   /* up to three genotypes possible */
   for (i=0; i<=2; i++)
 	 if (i != gmiss) for (j=0; j<=2; j++)
       if (j != gmiss) 
	   {
		  pr = log_condProb[1-i][1-j];

		  for (k=0, denom=0.0; k<=2; k++) if (k!=gmiss) 		   
		  {
             p = transProbL[k][i] * transProbR[j][k];
		     pr[1-k] = p;
		     denom += p;   /* calculate normalizing constant */
		  }

		  /* now normalize */
		  for (k=0; k<=2; k++) 
 			 pr[1-k] = log(pr[1-k]/denom);
	   }
}













/*===============================================================================

  Modified version of subroutines from MissMark.c (part of QTL Cartographer)

  ===============================================================================*/

/*
Assume that r and  a  are 3x3 square matrices.
Put a into r
*/
void AinR(double (*rr)[3], double (*aa)[3])
{
  int ii,jj;
  for ( ii = 0 ; ii<=2 ;ii++ )
    for ( jj = 0 ; jj<=2 ; jj++ )
      rr[ii][jj] = aa[ii][jj];
}  

/*
Assume that r, a, and b are 3x3 square matrices.
Calculate r = a.b
*/
void AdotB(double (*rr)[3], double (*aa)[3], double (*bb)[3])
{
  int ii,jj,kk;
  for ( ii =0  ; ii<=2  ;ii++ )
    for ( jj = 0 ; jj<=2  ; jj++ ) 
      for ( kk = 0  ; kk <= 2 ; kk++ )
        rr[ii][jj] = rr[ii][jj] + aa[ii][kk]*bb[kk][jj];
}  




/*
  Calc. conditional probability vector
  
  Return the conditional probability vector
  of the AA, Aa and aa genotypes in quk.
*/

void cond_prob(params *theparams, double (*Muv)[3], double *quk, double *pRuk, double *pLuk)
{
  double pCp[3],denom;

  int i;
  
  for ( i=0;i<3;i++ )
    pCp[i] = pRuk[i]*pLuk[i];
  for ( i=0;i<3;i++ )
    pRuk[i] = quk[i]*pCp[i];
  denom = 0;
  for ( i=0;i<3;i++ )
    denom = denom + quk[i]*pCp[i];
  if (denom == 0.0 ) {
      printf("\ndenom = 0.0 in cond_prob...");
    denom = 1.0;
  }
  for ( i=0;i<3;i++ )
    pLuk[i] = pRuk[i]/denom;
  a_Mv(quk,Muv,pLuk);

}



void markov_rec(params *theparams, double (*mm)[3], double rec, int rl)
{
  double mrec[3][3],tot[3][3],sftpm[3][3],trec,twrec;
  int i,ii,jj;
  for ( ii = 0 ; ii<3 ;ii++ )
    for ( jj = 0 ; jj<3 ; jj++ )
      tot[ii][jj] = mrec[ii][jj] = 0.0;
      
  if ( theparams->cross == 1 ) { /* B1 */
    mrec[1][1] = mrec[0][0] = 1.0 - rec;
    mrec[0][1] = mrec[1][0] = rec;    
  }
  else if ( theparams->cross == 2 ) { /* B2 */
    mrec[1][1] = mrec[2][2] = 1.0 - rec;
    mrec[2][1] = mrec[1][2] = rec;
  }
  else if ( theparams->cross == 3 ) { /* SFx or SFx crossed to P1 or P2 or Design III*/
    if ( theparams->tcross == 0 || theparams->tcross == 12 || theparams->tcrosst > theparams->crosst )
      selfed_f_tpm(mrec,rec,theparams->crosst);
    else if (theparams->tcross == 1) /* test cross to P1 */
      bselfed_f_tpm(mrec,rec,theparams->crosst,1);
    else if (theparams->tcross == 2) /* test cross to P2 */
      bselfed_f_tpm(mrec,rec,theparams->crosst,2);

/*  If there were more than one generations of backcrossing from the SFx, then
     we want to modify the G matrix for that number.    */
  if ( theparams->tcross < 3 && theparams->tcrosst > 1 ) { 
    if ( theparams->tcross == 1 ) { /* B1 */
      sftpm[1][1] = sftpm[0][0] = 1.0 - rec;
      sftpm[0][1] = sftpm[1][0] = rec;    
    }
    else if ( theparams->tcross == 2 ) { /* B2 */
      sftpm[1][1] = sftpm[2][2] = 1.0 - rec;
      sftpm[2][1] = sftpm[1][2] = rec;
    }
      AinR(tot,mrec);
      for ( i = 2 ; i<= theparams->crosst ; i++ ) {
        AdotB(mrec,tot,sftpm);
        AinR(tot,mrec);
      }
    for ( ii = 0 ; ii<3 ;ii++ )
      for ( jj = 0 ; jj<3 ; jj++ )
        tot[ii][jj] = sftpm[ii][jj] = 0.0;    
        
  }



  }
  else if ( theparams->cross == 4 ) { /* RFx */
    if ( theparams->tcross == 0 ) { /* no test cross */
      if ( theparams->crosst == 2 )
        trec = rec;
      else
        trec = 0.5 * (1.0 - pow( 1.0-rec, (double)(theparams->crosst-2) ) * (1.0-2.0*rec) );
      twrec = 1.0-trec;
      mrec[0][0] = mrec[2][2] = twrec*twrec;
      mrec[2][0] = mrec[0][2] = trec*trec;
      mrec[1][0] = mrec[1][2] = trec*twrec;
      mrec[0][1] = mrec[2][1] = 2.0*trec*twrec;
      mrec[1][1] = twrec*twrec + trec*trec;
    }
    else { /*   test cross to P1 or P2  */
      trec = 0.5 * (1.0 - pow( 1.0-rec, (double)(theparams->crosst-1) ) * (1.0-2.0*rec) );
      twrec = 1.0-trec;
      mrec[0][0] = mrec[1][1] = mrec[2][2] = 1.0-twrec ;
      mrec[2][0] = mrec[0][2] = 0.0;
      mrec[1][0] = mrec[1][2] = mrec[0][1] = mrec[2][1] = trec ;
      if ( theparams->tcross == 1 )
        mrec[2][1] = mrec[2][2] = mrec[1][2] = 0.0;
      else if ( theparams->tcross == 2 )
        mrec[0][0] = mrec[0][1] = mrec[1][0] = 0.0;
    }
  }
  else if ( theparams->cross == 5 ) {   /* RI lines...0=>Trudy, 1=>selfed RI, 2=>sibbed RI */
    mrec[0][0] = mrec[2][2] = 1.0 - rec;
    mrec[0][2] = mrec[2][0] = rec;
    if ( theparams->crosst == 1 ) {
      mrec[0][0] = mrec[2][2] = 1.0/(1.0+2.0*rec);
      mrec[0][2] = mrec[2][0] = 2.0*rec/(1.0+2.0*rec);
    }
    else if ( theparams->crosst == 2 ) {
      mrec[0][0] = mrec[2][2] = 1.0/(1.0+6.0*rec);
      mrec[0][2] = mrec[2][0] = 4.0*rec/(1.0+6.0*rec);
    }
  }
/*  If there were more than one generations of backcrossing from the F1, then
     we want to power up the G matrix for that number.    
      
             */
  if ( theparams->cross < 3 && theparams->crosst > 1 ) { 
      AinR(tot,mrec);
      AinR(sftpm,mrec);
      for ( i = 2 ; i<= theparams->crosst ; i++ ) {
        AdotB(mrec,tot,sftpm);
        AinR(sftpm,mrec);
      }
    for ( ii = 0 ; ii<3 ;ii++ )
      for ( jj = 0 ; jj<3 ; jj++ )
        tot[ii][jj] = sftpm[ii][jj] = 0.0;    
        
  }
  switch (rl) {
    default :
    case -3 : 
      break;
    case -2 : /*a- genotypes indicate that AA are impossible, so 0 first column. */
      for (i=0;i<3;i++) 
        mrec[i][0] = 0.0;
      break;
    case -1 : /*aa genotypes indicate that AA, Aa are impossible, so 0 first 2 columns. */
      for (i=0;i<3;i++) 
        mrec[i][0] = mrec[i][1] = 0.0;
      break;
    case 0  : /*Aa genotypes indicate that AA, aa are impossible, so 0 first and third columns. */
      for (i=0;i<3;i++) 
        mrec[i][0] = mrec[i][2] = 0.0;
      break;
    case 1  : /*AA genotypes indicate that aa, Aa are impossible, so 0 second 2 columns. */
      for (i=0;i<3;i++) 
        mrec[i][2] = mrec[i][1] = 0.0;
      break;
    case 2  : /*A- genotypes indicate that aa are impossible, so 0 last column. */
      for (i=0;i<3;i++) 
        mrec[i][2] = 0.0;
      break;
  }
  AdotB(tot,mm,mrec);
  AinR(mm,tot);
}




void calc_trans_prob(params* theparams, double rec,double** mrec)
{
  /* mrec[i][j] is prob of marker having geno j given a marker
	 'rec' (recombination rate, 0 to 0.5) away has genotype i
	 
	   0 => AA
	   1 => Aa
	   2 => aa
  */

  double tot[3][3],sftpm[3][3],trec,twrec;
  int i;
  int g0, g1, gmiss;
  

  if ( theparams->cross < 3)     /* B1 or B2 */
  {
      g1=1;
	  g0= (theparams->cross == 1)? 0:2;  /* B1 (cross=1) to parent 1 (AA) has g0=0 */
	  gmiss = 2-g0;                      /* the genotype not possible in this cross */
	                                        

	  if (theparams->crosst <= 1) 
	  {
		  mrec[g0][g0]=mrec[g1][g1] = 1-rec;
	      mrec[g0][g1]=mrec[g1][g0] = rec;
	  }
	  else 
	  {
		  trec = pow(1-rec,theparams->crosst);
		  mrec[g0][g0] = trec;
		  mrec[g0][g1] = 1 - trec;
		  mrec[g1][g0] = 2*pow(0.5,theparams->crosst)*mrec[g0][g1];
		  mrec[g1][g1] = 1 - mrec[g1][g0];
	  }
	  for (i=0; i<=2; i++) mrec[gmiss][i]=mrec[i][gmiss]=0.0;
  }
  else if (theparams->cross == 5)
  {
	  if (theparams->crosst==0)
		  mrec[0][2] = mrec[2][0] = rec;
	  else if (theparams->crosst == 1 )
  	      mrec[0][2] = mrec[2][0] =  2.0*rec/(1.0+2.0*rec);
	  else
          mrec[0][2] = mrec[2][0] =  6.0*rec/(1.0+6.0*rec);
          
	 mrec[0][0] = mrec[2][2] = 1 - mrec[0][2];
     for (i=0; i<=2; i++) mrec[1][i]=mrec[i][1]=0.0;
  }
  else if ( theparams->tcross == 0 && theparams->crosst == 2 )
  {
	  twrec = 1-rec;
      mrec[0][0] = mrec[2][2] = twrec*twrec;
      mrec[2][0] = mrec[0][2] = rec*rec;
      mrec[1][0] = mrec[1][2] = rec*twrec;
      mrec[0][1] = mrec[2][1] = 2.0*rec*twrec;
      mrec[1][1] = twrec*twrec + rec*rec;
  }
  else if ( theparams->cross == 3 ) { /* SFx or SFx crossed to P1 or P2 or Design III*/
    if ( theparams->tcross == 0 || theparams->tcross == 12 || theparams->tcrosst > theparams->crosst )
      selfed_f_tpm2(mrec,rec,theparams->crosst);
    else if (theparams->tcross == 1) /* test cross to P1 */
      bselfed_f_tpm2(mrec,rec,theparams->crosst,1);
    else if (theparams->tcross == 2) /* test cross to P2 */
      bselfed_f_tpm2(mrec,rec,theparams->crosst,2);

/*  If there were more than one generations of backcrossing from the SFx, then
     we want to modify the G matrix for that number.    */
    if ( theparams->tcross < 3 && theparams->tcrosst > 1 ) { 
      if ( theparams->tcross == 1 ) { /* B1 */
        sftpm[1][1] = sftpm[0][0] = 1.0 - rec;
        sftpm[0][1] = sftpm[1][0] = rec;    
	    for (i=0; i<=2; i++) sftpm[i][2]=sftpm[2][i]=0;
	  }
      else if ( theparams->tcross == 2 ) { /* B2 */
        sftpm[1][1] = sftpm[2][2] = 1.0 - rec;
        sftpm[2][1] = sftpm[1][2] = rec;
	    for (i=0; i<=2; i++) sftpm[i][0]=sftpm[0][i]=0;
	  }
  
      for ( i = 2 ; i<= theparams->tcrosst ; i++ ) {
        AinR2(tot,mrec);
        AdotB2(mrec,tot,sftpm);
      }    
	}
  }
  else if ( theparams->cross == 4 ) { /* RFx */
    if ( theparams->tcross == 0 ) { /* no test cross */
      trec = 0.5 * (1.0 - pow( 1.0-rec, (double)(theparams->crosst-2) ) * (1.0-2.0*rec) );
      twrec = 1.0-trec;
      mrec[0][0] = mrec[2][2] = twrec*twrec;
      mrec[2][0] = mrec[0][2] = trec*trec;
      mrec[1][0] = mrec[1][2] = trec*twrec;
      mrec[0][1] = mrec[2][1] = 2.0*trec*twrec;
      mrec[1][1] = twrec*twrec + trec*trec;
    }
    else { /*   test cross to P1 or P2  */
      trec = 0.5 * (1.0 - pow( 1.0-rec, (double)(theparams->crosst-1) ) * (1.0-2.0*rec) );
      twrec = 1.0-trec;
      mrec[0][0] = mrec[1][1] = mrec[2][2] = 1.0-twrec ;
      mrec[2][0] = mrec[0][2] = 0.0;
      mrec[1][0] = mrec[1][2] = mrec[0][1] = mrec[2][1] = trec ;
      if ( theparams->tcross == 1 )
        mrec[2][1] = mrec[2][2] = mrec[1][2] = 0.0;
      else if ( theparams->tcross == 2 )
        mrec[0][0] = mrec[0][1] = mrec[1][0] = 0.0;
    }
  }
}


  
void genProbs(int nn, params* theparams, int bc, int gmiss, 
			  QTL_INFO* aqtl, int revjump) 
{
	int i, q;
	int* prevQ, *nextQ, *geno;
	double prob[3];
	double* log_prob;
	int prevRenew, nextRenew;
	mygenome* qptr = aqtl->qptr;
	int disableOptimize = ((revjump & MCMC_SELECT_OPTIMIZE) == MCMC_OPTIMIZE_NONE);
	int disableProb = (revjump & MCMC_OPTIMIZE_NONE);

	if (aqtl->probValid) return;

    
    if (setNewFlank(qptr, &prevRenew, &nextRenew,
		            &aqtl->prevDist, &aqtl->nextDist, revjump))
          setupTable(gmiss,theparams, prevRenew, nextRenew, aqtl);
		  	         

    setUpGeno(qptr, &prevQ, &nextQ, &geno);

    /* generate probabilities for geno given flanking markers (old) */
    for(i=1;i<=nn;i++, prevQ++, nextQ++)
	{
	   if (*prevQ <-1 || *prevQ >1 || *nextQ<-1 || *nextQ > 1 || 
		    disableOptimize) 
	   {
		  log_prob = aqtl->log_prob[i] = aqtl->missing_prob[i];
	      calc_cond_prob2(theparams, bc, aqtl->qptr, i, 
			              &prob[2], &prob[1], &prob[0]);
		  for (q=-1; q<=1; q++) if (q != gmiss) log_prob[q] = log(prob[q+1]);
		     
	   }
	   else 
 	      aqtl->log_prob[i] = aqtl->log_condProb[*prevQ][*nextQ]; 
	}

	aqtl->probValid = !disableProb;  
	           /* we've just generated it, so probability array valid ...
	              unless we disable this feature */
}





void calc_cond_prob2(params *theparams,int bc, mygenome *gptr,int kk,
					double *pAA,double *pAa,double *paa)
{
  double pLuk[3],pRuk[3],Muv[3][3],quk[3],mml[3][3],mmr[3][3],rec;
  mygenome *tgptr;
  int ii,jj,k,lmark,rmark,mark,go_on;
  lmark=rmark=-3;
  for ( ii = 0 ; ii < 3 ; ii++ )
    for ( jj = 0 ; jj < 3 ; jj++ )
      Muv[ii][jj]=mml[ii][jj]=mmr[ii][jj] = 0.0;
  for ( ii = 0 ; ii < 3 ; ii++ )
    quk[ii]=Muv[ii][ii]=mml[ii][ii]=1.0;

/* For test crosses.  SFx then either Design III or a few more SF generations. */
  if ( theparams->cross == 3 && (theparams->tcross == 3 || theparams->tcross == 12 ) && theparams->tcrosst > theparams->crosst ) {
     Muv[1][1] = 1.0/pow(2.0, (double) (theparams->tcrosst - theparams->crosst));
     Muv[1][0] = Muv[1][2] = 0.5*(1.0-Muv[1][1]);
     Muv[0][0] = Muv[2][2] = 1.0;  /* This may make T(SFx)SFy work...but will it screw up D3's? */
     
    if (  theparams->tcross == 12 ) {
      if (  bc == 1 ) {
        mml[1][0] = mml[1][1] = 0.5;
        mml[2][2] = 0.0;
        mml[2][1] = 1.0;
      }
      else if ( bc == 2 ) {
        mml[1][2] = mml[1][1] = 0.5;
        mml[1][1] = 0.0;
        mml[0][1] = 1.0;
      }
      for ( ii = 0 ; ii < 3 ; ii++ )
        for ( jj = 0 ; jj < 3 ; jj++ ) 
          for ( k = 0 ;k<3;k++ )
            mmr[ii][jj] = mmr[ii][jj]+Muv[ii][k]*mml[jj][k];
      for ( ii = 0 ; ii < 3 ; ii++ )
        for ( jj = 0 ; jj < 3 ; jj++ ) 
          Muv[ii][jj] = mmr[ii][jj];
    }
    quk[1] = 1.0/pow(2.0, (double) (theparams->crosst - 1));
    quk[0]=quk[2] = 0.5*(1.0-quk[1]);
  }
  else
    for ( ii = 0 ; ii < 3 ; ii++ )
      quk[ii]=Muv[ii][ii] =1.0;
  
  for ( ii = 0 ; ii < 3 ; ii++ )
    for ( jj = 0 ; jj < 3 ; jj++ )
       mml[ii][jj]=mmr[ii][jj] = 0.0;
  for ( ii = 0 ; ii < 3 ; ii++ )
    mml[ii][ii]=mmr[ii][ii]=1.0;
  *pAA = *pAa = *paa = 0.0;
  if ( gptr->markr > 0 )
    mark =  gptr->genotype[kk];
  else
    mark = -3;
  switch (mark) {
    case -1 :
      *paa = 1.0;
      break;
    case 0 :
      *pAa = 1.0;
      break;
    case 1 :
      *pAA = 1.0;
      break;
    case -2 :
      mml[0][0] = 0.0;
    case 2 :
      if (mark == 2)
        mml[2][2] = 0.0;
    default :
      tgptr = gptr;
      if ( tgptr->prev != NULL ) {
        tgptr = tgptr->prev;
        go_on = 1;
      }
      else 
        go_on = 0;
      while ( go_on == 1 ) 
        if ( tgptr->chrom == gptr->chrom ) {
          rec = mapfunc( tgptr->dist , -1);
          lmark =  tgptr->genotype[kk];
  
          markov_rec( theparams, mml,  rec  , lmark  );
          if ( lmark > -2 && lmark < 2 )
            go_on = 0;
          if ( tgptr->prev != NULL   ) 
            tgptr = tgptr->prev;
          else 
            go_on = 0;
        }
        else
          go_on = 0;
      a_Mv(pLuk,mml,quk);
            
      tgptr = gptr;
      if ( tgptr->next != NULL ) 
        go_on = 1;
       else 
        go_on = 0;
      while ( go_on == 1 ) 
        if ( tgptr->next->chrom == gptr->chrom ) {
          rec = mapfunc( tgptr->dist , -1);
          rmark =  tgptr->next->genotype[kk];

          markov_rec(theparams,mmr,  rec  , rmark  );
          if ( rmark > -2 && rmark < 2 )
            go_on = 0;
          tgptr = tgptr->next;
          if ( tgptr->next == NULL   ) 
            go_on = 0;
        }
        else
          go_on = 0;
      a_Mv(pRuk,mmr,quk); /* pRuk =  Hzk+1.Hzk+2...Hzl.c */
      assign_q(quk,theparams);
      cond_prob(theparams,Muv,quk,pRuk,pLuk);
      *pAA = quk[0];
      *pAa = quk[1];
      *paa = quk[2];    
      break;
  }
}




/*
Simple subroutine to do a = M.v, where a, v are 0 offset 3 vectors
and M is a 0,0 offset 3x3 matrix.
*/
void a_Mv(double *a, double (*M)[3], double *v)
{
int i,j;
  for ( i=0;i<3;i++ ) {
    a[i] = 0.0;
    for ( j=0;j<3;j++ )
      a[i] = a[i] + M[i][j]*v[j];
  }
}





/*

This assigns the vector qk' which is written explicitly for an F2
population in equation 3 of Jiang and Zeng's paper on dominanant and
missing markers.

*/
void assign_q(double *ql,params *theparams)
{
  int ii;
  for ( ii = 0 ; ii < 3 ; ii++ )
    ql[ii] = 0.0;
  if ( theparams->cross == 1 ) {  /*  pr(AA) = 1-pr(Aa);  pr(Aa) = pow(.5,t) for B1t */
    ql[1] = pow(0.5, (double) theparams->crosst) ;
    ql[0] = 1.0 - ql[1];
  }
  else if ( theparams->cross == 2 ) {  /*  pr(aa) = 1-pr(Aa);  pr(Aa) = pow(.5,t) for B2t */
    ql[1] = pow(0.5, (double) theparams->crosst) ;
    ql[2] = 1.0 - ql[1];
  }
  else if ( theparams->cross == 3 ) {  /*  pr(AA) = pr(aa) = .25 for B1 */
    if ( theparams->tcross == 1 )
      ql[0]=ql[1]=0.5;
    else if  ( theparams->tcross == 2 )
      ql[2]=ql[1]=0.5;
    else  {
      ql[1] = 1.0/pow(2.0, (double) (theparams->crosst-1));
      ql[0] = ql[2] = 0.5*(1.0-ql[1]);
    } 
  }
  else if ( theparams->cross == 4 ) {
    if ( theparams->tcross == 1 )
      ql[0]=ql[1]=0.5;
    else if  ( theparams->tcross == 2 )
      ql[2]=ql[1]=0.5;
    else  {
      ql[0] = ql[2] = 0.25 ;
      ql[1] = 0.5;
    }  
  }
  else if ( theparams->cross == 5 ) 
    ql[0] = ql[2] = 0.5 ;
}

/*
  Suppose we are producing Selfed Ft's.  This will calculate the
  transition matrix.
*/
void selfed_f_tpm(double (*sftpm)[3], double r, int t)
{
  double twotm2,twotm1, c1, c2, c3;

  if ( t != 2 )
    twotm2 = pow( 2.0, (double) t - 2.0 );
  else
    twotm2 = 1.0;
  twotm1 = twotm2*2.0;
  c1 = 1.0+2*r;
  c2 = pow( (0.5-r), (double) t ) ;
  c3 = pow( 0.5-r*(1.0-r) , (double) t - 1.0 );
  
  sftpm[0][0] = sftpm[2][2] =   (twotm1/c1       - 1.0 - twotm1*c2/c1 + twotm2*c3)/(twotm1-1.0);
  sftpm[0][1] = sftpm[2][1] =   (                  1.0 -                twotm1*c3)/(twotm1-1.0);
  sftpm[0][2] = sftpm[2][0] =   (2.0*r*twotm1/c1 - 1.0 + twotm1*c2/c1 + twotm2*c3)/(twotm1-1.0);
  sftpm[1][0] = sftpm[1][2] =   0.5 -  twotm2*c3 ;
  sftpm[1][1] =                 twotm1*c3 ;

}

void bselfed_f_tpm(double (*sftpm)[3], double r, int t, int bc)
{
  double c1, c2;


  c1 = pow( (0.5-r), (double) t);
  c2 = 1.0+2.0*r;
   
  sftpm[0][0] = sftpm[1][1] = sftpm[2][2] =  c1 + (1.0-c1)/c2;
  sftpm[0][1] = sftpm[1][0] = sftpm[0][2] = sftpm[2][1] =  (2.0*r+c1)/c1 - c1;
  if (bc==1) 
    sftpm[0][2] = sftpm[1][2] = sftpm[2][2] = sftpm[2][0] = sftpm[2][1] =  0.0;
  if (bc==2)
    sftpm[0][0] = sftpm[0][1] = sftpm[0][2] = sftpm[1][0] = sftpm[2][0] =  0.0;
}



void selfed_f_tpm2(double** sftpm, double r, int t)
{
  double twotm2,twotm1, c1, c2, c3;

  if ( t != 2 )
    twotm2 = pow( 2.0, (double) t - 2.0 );
  else
    twotm2 = 1.0;
  twotm1 = twotm2*2.0;
  c1 = 1.0+2*r;
  c2 = pow( (0.5-r), (double) t ) ;
  c3 = pow( 0.5-r*(1.0-r) , (double) t - 1.0 );
  
  sftpm[0][0] = sftpm[2][2] =   (twotm1/c1       - 1.0 - twotm1*c2/c1 + twotm2*c3)/(twotm1-1.0);
  sftpm[0][1] = sftpm[2][1] =   (                  1.0 -                twotm1*c3)/(twotm1-1.0);
  sftpm[0][2] = sftpm[2][0] =   (2.0*r*twotm1/c1 - 1.0 + twotm1*c2/c1 + twotm2*c3)/(twotm1-1.0);
  sftpm[1][0] = sftpm[1][2] =   0.5 -  twotm2*c3 ;
  sftpm[1][1] =                 twotm1*c3 ;

}


void bselfed_f_tpm2(double **sftpm, double r, int t, int bc)
{
  double  c1, c2;


  c1 = pow( (0.5-r), (double) t);
  c2 = 1.0+2.0*r;
   
  sftpm[0][0] = sftpm[1][1] = sftpm[2][2] =  c1 + (1.0-c1)/c2;
  sftpm[0][1] = sftpm[1][0] = sftpm[0][2] = sftpm[2][1] =  (2.0*r+c1)/c1 - c1;
  if (bc==1) 
    sftpm[0][2] = sftpm[1][2] = sftpm[2][2] = sftpm[2][0] = sftpm[2][1] =  0.0;
  if (bc==2)
    sftpm[0][0] = sftpm[0][1] = sftpm[0][2] = sftpm[1][0] = sftpm[2][0] =  0.0;
}




/*
Assume that r and  a  are 3x3 square matrices.
Put a into r
*/
void AinR2(double (*rr)[3], double **aa)
{
  int ii,jj;
  for ( ii = 0 ; ii<=2 ;ii++ )
    for ( jj = 0 ; jj<=2 ; jj++ )
      rr[ii][jj] = aa[ii][jj];
}  


/*
Assume that r, a, and b are 3x3 square matrices.
Calculate r = a.b
*/
void AdotB2(double **rr, double (*aa)[3], double (*bb)[3])
{
  int ii,jj,kk;
  for ( ii =0  ; ii<=2  ;ii++ )
    for ( jj = 0 ; jj<=2  ; jj++ ) 
      for ( kk = 0  ; kk <= 2 ; kk++ )
        rr[ii][jj] = rr[ii][jj] + aa[ii][kk]*bb[kk][jj];
}  
