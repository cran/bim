#ifndef REVJUMP_GUARD
#define REVJUMP_GUARD
/**************************************************************
  File:         revjump.h
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:
  -------

    Header file for MCMC program.  See also chrom.h and qtl.h

**************************************************************/


#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>


#include "chrom.h"
#include "qtl.h"
#include "mygenome.h"
#include "Utilities.h"
#include "R.h"




#include <assert.h>



typedef struct Params {  /*Structure to hold parameters*/
  int cross;       /*primary cross type */
  int crosst;      /*primary cross generations   */
  int tcross;      /*type of test cross*/
  int tcrosst;     /*generations of test cross*/
  int chrom;
  int seed;
  char* error;
} params;



typedef struct Data 
{
  /* data for analysis */
  int nn;
  int nChrom;
  int totalMark;
  double totalChromLen;
  double* y;
  int bc;
  double ybar;
  double y_var;
  double y_stdev;
  params* theparams;

  CHROMOSOME* myChroms;
  double* chrom_pos;
   

  /* death-birth stuff */
  double* bp;
  double* dp;
  double* prior_ratio;

  /* some parameters */
  int nQtl;
  double mu;
  double sigmasq;

  /* cross parameters */
  int gmiss;
  int na[3];

  QTL_INFO** myQtls;

} DATA;



int qtl_compare(const void* a, const void* b);

#define MAX_INTERACTION 100
#ifndef MAXQTL
#define MAXQTL 30
#define MAX_CHOL (MAXQTL*2+1)
#endif





/* used to index prior mean,var and log_var (mean, additive and dominance respectively */
#define MU  QTL_NONE
#define ADD  QTL_ADD
#define DOM  QTL_DOM 


#define MAX_PARAM 2
#define MAX_ADD_PARAM 2



#define INFINITY 1E40

/* codes for birth/death moves */
#define C_SWAP 5
#define C_EXCHANGE 4
#define C_UPDATE 3
#define C_BIRTH 1
#define C_DEATH 2

#define NOT_IMPLEMENTED 4
#define BAD_CODE 5




/*--------------------------------------------------------------------------
 * One field in the MCMC_PARAM field contains a code which controls the 
 * analysis prefereence.  Generally it is set to 1, which activates RJ-MCMC,
 * a Poisson prior on the number of QTL, JAYA sampling for new effects (in 
 * RJ-MCMC).  Also additive/dominance effects are sampled as dictated by the
 * design.
 *--------------------------------------------------------------------------
 */

#define MCMC_SELECT_PRIOR_FIELD   0xF      
/* distrib for prior number of QTL  */
#define MCMC_SELECT_SAMPLE_FIELD  0xF0     
/* method in RJ-MCMC of proposing effect */
#define MCMC_SELECT_DOM_FIELD     0xF00    
/* obtain estimates for dominance        */
#define MCMC_SELECT_EFFECT_FIELD  0xF000   
/* type of variance on prior for effects */
#define MCMC_SELECT_PROPOSE_FIELD 0xF0000  
/* proposal from birth death is drawn 
                                    from this distribution*/
#define MCMC_FLAGS                0xF00000 
/* bit flags, described below, to turn
                                     off default features */
#define MCMC_SELECT_MIX           0x3000000
/* change mix from 90-10 birth/death to
								     update (2 bits used), 3rd bit also used as
                                     a flag */
#define MCMC_SELECT_OPTIMIZE      0x30000000  
/* some optimization code bit flags (2 bits used) */
#define MCMC_SELECT_CHROM         0xC0000000                        
                                           


/***************************/
/* MCMC_SELECT_PRIOR_FIELD */
/***************************/
#define SELECT_SIMPLE_MCMC      0         
#define SELECT_NQTL_POISSON     1         
#define SELECT_NQTL_GEOMETRIC   2
#define SELECT_NQTL_UNIFORM     3    
#define SELECT_NQTL_FISCH       4

/****************************/
/* MCMC_SELECT_SAMPLE_FIELD */
/****************************/
#define SELECT_SAMPLE_DECOMP      0x00    /* default, use decomposition */
#define SELECT_SAMPLE_FISCH       0x10 

/**************************/
/* MCMC_SELECT_DOM_FIELD  */
/**************************/
#define SELECT_GET_APPROP         0x000    /* default,get appropriate model */
#define SELECT_GET_ADD_ONLY       0x100  
#define SELECT_GET_DOM_ONLY       0x200    
#define SELECT_RJMCMC             0x300   /* select components to add 
                                             using RJ-MCMC */

/****************************/
/* MCMC_SELECT_EFFECT_FIELD */
/****************************/
#define SELECT_VAR_GAFFNEY        0x0000
#define SELECT_VAR_FISCH          0x1000  

/****************************/
/* MCMC_SELECT_PROPOSE      */
/****************************/
#define SELECT_PROPOSE_PRIOR      0x00      /* default, propose birth-death 
                                               probs from prior */
#define SELECT_PROPOSE_POISSON    0x10000   /* always propose poisson */
#define SELECT_PROPOSE_GEOMETRIC  0x20000   /* always propose geometric */
#define SELECT_PROPOSE_UNIFORM    0x30000   /* always propose uniform */
#define SELECT_PROPOSE_FISCH      0x40000   /* always propose geometric */


/***************************/
/* MCMC_FLAGS              */
/***************************/
/* these flags allow additional options to be set or not (one bit).
   The previous ones allow a range of options (up to 16, x0 to xF,
   4 bits worth)  */

#define MCMC_SELECT_100PCNT_UPDATE 0x100000  /* default update parameters
                                                only if no birth or 
                                                death update */
                                                
#define MCMC_SELECT_BD_NOTEXACT    0x200000  /* here bp+bp <= cval.
                                                Whereas, by default 
                                                we have bp+cp=cval.
                                                cval is set by setting one of
                                                the MCMC_OPTIMIZE flag in 
                                                'revjump' */   
											
#define MCMC_SELECT_SIMPLE_GIBBS   0x400000  /* default is to use a multiple 
                                                parameter update for effects 
                                                (rather than simple Gibbs) */
											          
#define MCMC_SELECT_MOVE_NONE      0x800000    /* default is to propose 
                                                  update on any chrom */

/***************************/
/* MCMC_MIX                */
/***************************/
/* 2 bit field controlling birth-death-update mix */
#define MCMC_OPTIMIZE_DEFAULT     0x0000000    /* standard 90-10 */
#define MCMC_OPTIMIZE_80_20       0x1000000    /* standard 80-20 */
#define MCMC_OPTIMIZE_70_30       0x2000000    /* standard 70-30 */
#define MCMC_OPTIMIZE_50_50       0x3000000    /* standard 50-50 */


/* set this flag to always turn on long-range update (regardless of
  whether a birth/death occured.  Default is to have a long-range
  jump only when updating QTL genotypes */
#define MCMC_OPTIMIZE_BOOST       0x4000000    


/***************************/
/* MCMC_OPTIMIZE           */
/***************************/
/* 2 bit field controlling optimization -- 3 bits*/
#define MCMC_OPTIMIZE_NOPROB     0x10000000  /* turn off prob */
#define MCMC_OPTIMIZE_NOTABLE    0x20000000  /* turn off prob andtable  */
#define MCMC_OPTIMIZE_NONE       0x30000000  /* generate each individual's
                                                prob of QTL geno given
                                                marker */
#define MCMC_FLAG_RANDOM_CHROM   0x40000000  /* default is to select the 
						chromosome with prob proportional
						to its length.  This flag makes
						each chromosome equally likely
						to be chosen (affects acceptance
						ratio) ... lowers it for smaller
						chromosomes.
												*/
#define MCMC_FLAG_NEW_PRIOR   0x50000000  /* selects each chromosome with equal
					     probability, and then position is 
					     uniform over length.  However, this
					     is considered the proior (rather than
					     uniform over the genome).
					  */








typedef struct workData
{
  /* effect and proposed (new) effect */

  /* used in QR decomposition of matrix of genotypes X
	 * with Y = phenotype values, Z = proposed (new) geno,
	 *
	 *  K = Z'RY,   M= 1/(Z'RZ), QZ = Q'Z, QY = Q'Y
	 *  R = I - X{(X'X)^-1}X' = I - QQ', 
	 *
	 *  Note: K = (QZ)'(QY) and M = 1/(QZ'QZ)
	 *        See get_effect(...) in birth.c for full discussion
     */	
  double* u;        /* proposed uniform number */
  double* resid;    /* residual (y-mu-effect) for all individuals */
  double* newResid; /*?????*/
  int* perm_num;   /* use to permute QTL in random_perm_QTL */
  
  /* used in update_params.c */
  double* XtY;
  double* new_XtY;
  double** XtX;
  double** new_XtX;
  double** chol;
  double* p;
  double* pvar;
  double* pmean;
  double* mod_effect;
  double* weight;
  double* work;
  double* oldWeight;   /* needed in proposeFischDeath */
  double* oldVar;      /* needed in proposeFischDeath */

  /* all we need to retain from effect update (for birth/death acceptance ratio */
  double bCb;
  double d_invD_d;
  double log_det_Chol;
  double log_det_invD;


  double* new_r;
  double* log_norm_const;

  QTL_INFO** output_qtls;

  double* minProbs;         /* for diagnostics */
  double* avgMinProbs; 



} WORK;


typedef struct MCMCParam
{
  /* 1 => reversible jump, 0 otherwise. */
  int revjump;
  int niter;
  int nby;
  double cval;

  int offset;
  /*  int numParam;*/    /* number of parameters to add in birth/death step */
  int addParam;

  double burnIn;
  double preBurnIn;

  double HM;
  double SHM;
  int idx;

} MCMC_PARAM;




typedef struct priorstruct
{
  double mean[3];
  double var[3];
  double alpha[3];          /* prior for variance for add/dom var, 
			       which is Beta(1,beta) */
  double beta[3];          /* prior for variance for add/dom var, 
			      which is Beta(1,beta) */
  double log_var[3];
  int sampleVar[3];


  double log_mu_var;
  double sig_a1;
  double sig_a2;
  double qtl_mean;
  double p;

  int priorDistribution;          /* can be one of  POISSON   Here qtl_mean is the mean of the distrib.UNIFORM   Here qtl_mean holds the maximum number of QTL GEOMETRIC  "   qtl_mean is the mean number of QTL (1/q) Prob of x QTLs = q^x * p (x=0,1,2,...) where q=1-p MIXED     Here qtl_mean is the cutoff before which the prior is uniform.  After there is an geometic decay at rate p. */
} PRIORS;




typedef struct Missing {
  int lMark;     /* 0 => undefined, otherwise there is a valid flanking marker */        
  int rMark;     

  double lPos;   /* position of marker */
  double rPos;

  double scale;  /* In general, we can compute probabilities of a qtl [given
		    flanking markers] to within a constant of proportionality.
		    Only in simple cases (no A- or a- genotypes) can we
		    determine this constant uniquely */

  double** lProb;  /* for each missing marker, we have a vector of 3 values */    
  double** rProb;  /* this give prob of each genotype at the missing marker */    
} MISSING;


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>




#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif


#ifndef MIN
#define MIN(A, B) ((A) > (B) ? (B) : (A))
#endif




#ifndef ABS
#define ABS(A) ( ((A)>0)? (A) : -(A) ) 
#endif

#ifndef EPSILON
#define EPSILON 1E-7
#endif


/* Jaya ususally gives A=genunf() and B=prob[0], then
   returns a marker of genotype 1 wp B and -1 wp (1-B) */
#define random_marker(A,B) ((A) < (B) ? (1) : (-1))
#define ismiss(A) (((A) < -2) ? (1) : (0))


/* mcmc.c */
void mcmc(MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors,
          CHROMOSOME* chrInfo, WORK* myWork, double* result);
void initVars(int nn, double** nmove, int* num_qtl, MCMC_PARAM* myMCMC,
	      double qtl_prior_mn, double cval);
void setupBurnIn (MCMC_PARAM* myMCMC, DATA* myData);
void birth_death(int move_type, DATA* myData, MCMC_PARAM* myMCMC, 
		 WORK* myWork, PRIORS* priors, 
		 CHROMOSOME* chrInfo, double** nmove);
void update_effect(int nn, int nqtl, QTL_INFO* aqtl, QTL_INFO** all_qtls,
		   MCMC_PARAM* myMCMC, DATA* myData, 
		   PRIORS* priors, WORK* myWork);
void update_effects(int nn, int nQtl, int revjump, double* y, 
		    double* mu, double sigmasq, int* na, 
		    QTL_INFO** qtls, PRIORS* priors,WORK* myWork);
					  
					  

void normalProb(int nn, double* resid, double sigmaSq, 
		double* diag_fY, double* maxResid, double* pcnt_gt_99);
 void outputResults(int nn,
		   int iter, DATA* myData, QTL_INFO** all_qtls, 
		   double* resid, QTL_INFO** output_qtls, 
		   PRIORS* priors, MCMC_PARAM* myMCMC,
		   double* result, int* output_idx);


void outputStatistics1(int nn, MCMC_PARAM* myMCMC, params* theparams, 
		       PRIORS* priors, DATA* myData);
void outputStatistics2(int nn, MCMC_PARAM* myMCMC, params* theparams, 
		       double** nmove, int* num_qtl,
		       int countedIter, PRIORS* priors, DATA* myData);
double calc_h2(int nn, double* y, double mu, double sigmasq, double* resid);
void diagnose(int nn, QTL_INFO** all_qtls, WORK* myWork, 
	      DATA* myData, MCMC_PARAM* myMCMC);
void calcMeanVar(int nn, double* vals, double* mean, double* var);
void updateMean(int i, double val, double* mean);


/* lod.c */
double lodnull(int nn, double y_var);
double get_lod(int nn, double sigmasq, double ybar, double y_var,
	       double* resid);
/**********************************************************************
    **********************************************************************/



/* birth.c */


void setWeights(int nQtl, int* na,
		int revjump, QTL_INFO** all_qtls, 
		double* w, QTL_INFO* newQtl);

void get_new_locus(int revjump, 
		   int nChrom, double* chrom_pos, double totalChromLen, 
		   CHROMOSOME* chromInfo, CHROMOSOME** chrom, 
		   int* lmark, double* pos);
void get_local_locus(CHROMOSOME* chrom, int* lmark, double* pos);

void get_new_qtl_genotype(int nn, QTL_INFO* qtl);

double get_effect(long int nn, int nQtl, 
                  double* y, double mu, double sigmasq,
		  QTL_INFO* oldQtl, QTL_INFO* newQtl, 
		  QTL_INFO** qtls, int revjump,
		  PRIORS* priors, 
		  int* old_na, double* mod_effect,
		  double** XtX, double* XtY, double** chol, double* p, 
		  double* pmean, double* pvar, double* u, double* weight,
		  double* bCb, double* d_invD_d, 
                  double* log_det_Chol, double* log_det_invD,
		  WORK* myWork);
				           


QTL_INFO* createQtl(int nn, int qtlNum, QTL_INFO** p_newQtl, 
		    CHROMOSOME* chrom, int lmark, double new_position,
		    int addParam, double* a, double* w);
void setEffect(int nn, int nQtl, double* y, double* mod_effect, 
	       QTL_INFO** all_qtls, double* mu, 
	       double* w, double* resid, int* na);
void setCholParams(WORK* myWork, double bCb, double d_invD_d, 
		   double log_det_Chol, double log_det_invD);
void addQtlToChrom(QTL_INFO* aqtl);
double calcResidSS(int nn, double* resid);
double proposeFischBirth (long int nn, int nQtl, 
			  double* y, double mu, double sigmasq,
			  QTL_INFO* newQtl, 
			  QTL_INFO** qtls, int revjump, 
			  PRIORS* priors, int* na,
			  double* mod_effect,
			  double* pmean, double* pvar, double* weight,
			  double* oldVar, double* oldWeight,
			  double* resid, double* newResid);
double proposeFischDeath (long int nn, int nQtl, 
			  double* y, double mu, double sigmasq,
			  QTL_INFO** qtls, int revjump, 
			  PRIORS* priors, int* na,
			  double* mod_effect,
			  double* pmean, double* pvar, double* weight,
			  double* oldVar, double* oldWeight,
			  double* resid, double* newResid);
int getFischEffect(int nQtl, QTL_INFO** qtls, double mu,
		   PRIORS* priors, double* mod_effect, double* w, double* var);
QTL_INFO* swap_add_dom(int nn, int nQtl, QTL_INFO** qtls,
		       PRIORS* priors, 		              
		       DATA* myData, WORK* myWork, MCMC_PARAM* myMCMC);




/* death.c */
int death(DATA* myData, PRIORS* priors, 
	  params* theparams, WORK* myWork, MCMC_PARAM* myMCMC);
int get_death_interval(int nQtl);

void dropQtl(int nn, int* nQtl, double* y, QTL_INFO** all_qtls, 
	     double *mod_effect, double* w, double* mu, 
	     double* resid, int* na);
int removeQtlFromList(int nQtl, QTL_INFO* aqtl, QTL_INFO** qtls);
void moveQtlToEndofList(int nQtl, QTL_INFO* aqtl, QTL_INFO** qtls);
QTL_INFO* removeIQtlFromList(int nQtl, int i, QTL_INFO** qtls);


/* accept_birth.c */
int select_move(int nQtl, double* bp, double* dp);

void calcResid2(int nn, int nQtl, double* y, double* modified_parms,
		QTL_INFO** qtls, double* newResid);
double get_log_proposal_ratio(int nQtl, double* bp, double* dp, 
			      double* priorRatio);
double get_log_position_ratio(int revjump, CHROMOSOME* chrom, 
			      DATA* myData);
double get_prior_position_ratio(double totalChromLen, int nChrom, 
				int nMark, double IMlen);



/*  update_parms.c */
double update_mu(int nn, double mu, double sigmasq, double* resid,
		 double mu_prior_mn, double mu_prior_var);
double update_add_effect(int nn, int nQtl, double sigmasq,
			 QTL_INFO* modQtl, double* resid,
			 double effect_prior_mn, double effect_prior_var);
double update_dom_effect(int nn, int nQtl, double sigmasq,
			 QTL_INFO* modQtl, double* resid,
			 double effect_prior_mn, double effect_prior_var);
double Gibbs_update_sigmasq(int nn, double* resid,
			    double prior_sigsq_a1,
			    double prior_sigsq_a2);
double MH_update_sigmasq(int nn, double* resid,
			 double sigmasq,
			 double deltaS2, double maxS2);
double qtl_prior_effect_ratio (double effect_prop, double w_prop,
			       double effect, double w,
			       double effect_mean, double effect_var, 
			       double log_effect_var);
double calc_prior_effect_term(double effect, double w, double effect_mean,
			      double log_effect_var, double effect_var);

double update_effect_prior_var(int type, int nQtl,  QTL_INFO** all_qtls,
			       PRIORS* priors, double y_var);

void setUpGeno(mygenome* qptr, int** prevGeno, int** nextGeno, int** geno);
void swapRowAndCol(int nn, double** a, int row1, int row2);
void moveQtlToEndOfXtX(int nQtl, QTL_INFO** qtls, int changeQtl, 
		       double** XtX, double* XtY, int nterm);
void removeRowCol(int n, double** a, int row);
double log_determinant(int nterm, double* p, double sigmasq);
QTL_INFO* long_range_update(int nn, int nQtl, QTL_INFO** qtls, 
			    double* chrom_pos, double totalChromLen,
			    int nChrom, CHROMOSOME* chromInfo, PRIORS* priors, 		              
			    DATA* myData, WORK* myWork, MCMC_PARAM* myMCMC);



int cholesky(int n, double **a, double* p);
int incremental_cholesky(int n, double** chol, double* p, double * newcol);
void choleskySolve(int n, double **a, double* p, double* b, double* x);
void forwardSubs(int n, double **a, double*p, double* b);
void backwardSubs(int n, double **a, double*p, double* b);
void mycopy(int n, int m, double** to, double** from);
int setAddDomCovMatrix(int nn, int nQtl, QTL_INFO** selectedQtls, 
		       double* y, double** XtX, double* XtY, int* na);
				       
void addColToAddDom(int nn, int old_nQtl, int* old_na, 
		    QTL_INFO** qtls, QTL_INFO* aqtl,
		    double* y, double* XtY, double** XtX);

void setAddDomDiag_Row1(int nn, int idx, QTL_INFO* aqtl, 
			double* y, double** a, double* wr);
int generate_effects(int nterm, double** XtX, double* XtY, 
		     double* pmean, double* pvar,double sigma2, 
		     double** chol, double* p, 
		     double* u, double* mod_effect,
		     double* bCb, double* d_invD_d,
		     double* log_det_Chol, double* log_det_invD, int genBeta);
void calc_interaction(int idx, int idx2, double** a, 
		      QTL_INFO* aqtl, QTL_INFO* bqtl, int nn);
void setPriorMeanVar(int nQtl, int revjump, 
		     QTL_INFO** qtls, QTL_INFO* newQtl, PRIORS* priors, 
		     double* w, int* na, double* k, double* v);



/* update_qtl.c */
void setValidFlag(QTL_INFO* aqtl, int revjump);
double update_lambda_qtl(QTL_INFO* aqtl, QTL_INFO** all_qtls, 
			 DATA* myData, WORK* myWork, int revjump);
double get_new_lambda(int nMark, double lambda);
int setNewFlank(mygenome* qptr, int* prevRenew, int* nextRenew,
		double* prevDist, double* nextDist, int revjump);

void SwapGenome(mygenome** g1, mygenome** g2);
void Swap3DTable(double**** g1, double**** g2);
void Swap2DTable(double*** g1, double*** g2);
void SwapIVec(int** g1, int** g2);
void SwapDVec(double** g1, double** g2);
void SwapInt(int* g1, int* g2);
void SwapDble(double* g1, double* g2);
void swapQtl(QTL_INFO** g1, QTL_INFO** g2);
void swapQtlData(QTL_INFO* Qtl1, QTL_INFO* Qtl2);





double calc_flank_dist(mygenome* qptr);
void get_new_pos(CHROMOSOME* chrom, double oldPos, double delta, 
		 int* newMark, double* newPos);
double gen_qprob (double* a, int geno, double sigmasq,             
		  double resid,
		  double * prob, double *qprob, double* new_r, int gmiss);
void setupTable(int gmiss, params* theparams, 
		int renewPrev, int renewNext, QTL_INFO* aqtl);
void calc_trans_prob(params* theparams, double rec,double **mrec);
						
void normalizeProb(int gmiss, double*** condProb);
void setLambda(QTL_INFO* aqtl, int mark, double pos, int qtlNum);
void copySwapGenoTo(QTL_INFO* newQtl, QTL_INFO* aqtl);




/* initvals.c */
void GenGenotype(int gmiss, double* pr, int* genotype);
double IMdist(QTL_INFO* aqtl);
int* genotype1(QTL_INFO* aqtl);
int* igenotype(QTL_INFO* aqtl);
int** p_igenotype(QTL_INFO* aqtl);
void dgenotype(QTL_INFO* aqtl, int nn, double* dz);
int qtlnum(QTL_INFO* aqtl);
double haldane( double dist_val );
double gammln2(double xx);



/* read_data.c */
void getdata( MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, 
              CHROMOSOME* chromInfo, WORK* myWork,
	      int* revjump, double* burn_in, double* preburn_in, 
	      int* niter, int* nby, double* mu, double* sigmasq,
	      double* meanmu, double* sdmu, double* siga1, double* siga2,
	      double* meanadd, double *sdadd,double* meandom, double*sddom,
	      double* qtlmean);





/* qtl_move.c */
double calc_intermarker_dist(CHROMOSOME* chrom, int m);
mygenome* findInsertPos(double pos, mygenome* lmark, int chrom);
void insertQtl(mygenome* qptr, int lmark, double pos, 
	       CHROMOSOME* chrom, int qtlNum);
void removeQtl(mygenome* qptr);
int checkIntegrity(int nQtl, CHROMOSOME* chrom);
void replaceQtl(mygenome* old, mygenome* new);
void restoreQtl(mygenome* qptr);
int EQUALS(double x, double y);
int LE(double x, double y);
int LT(double x, double y);

/* numcmp.c */
int binSearch(int nval, double* vals, double searchVal);
int numcmp(double *v1, double *v2);



void checkResid(int nn, int nQtl, double mu, double* y, 
		QTL_INFO** allQtls, double* resid, int gmiss);

void calc_cond_prob2(params *theparams, int bc, mygenome *gptr,int kk,
		     double *pAA,double *pAa,double *paa);

void selfed_f_tpm2(double** sftpm, double r, int t);
void bselfed_f_tpm2(double **sftpm, double r, int t, int bc);
void AinR2(double (*rr)[3], double **aa);
void AdotB2(double **rr, double (*aa)[3], double (*bb)[3]);


void Cholesky(int n, double**z, double** L);
void Inverse(int n, double** z, double** inverse_z);
void Multiply(int n, int m,  double** a, double*b, double*c);
double Determinant(int n, double** z);




/* main.c */
void setupWork(int nn, WORK* myWork);
void setupBirthDeathProbs(int revjump, double qtl_param, 
			  double cval,
			  double** bp, double** dp, double** priorRatio);
void noramlizeBirthDeath(double* bp, double* dp, double cval);
void setupDiagnostics(int nn, MCMC_PARAM* myMCMC, params* theparams);



void random_perm_QTL(int nQtl, QTL_INFO** qtls, QTL_INFO** buff, 
		     int* perm_num);
void printX(int nn, int nQtl, QTL_INFO** qtls);
void printXtX(int nterm, double** XtX);





/* diagnostics */
void checkCholesky(int nQtl,int revjump, QTL_INFO** qtls, PRIORS* priors, 
		   double* weight, int* na, double* pmean,double* pvar,
		   double** XtX, double* XtY, double sigmasq,
		   double** chol, double* p, double* u, double* mod_effect,					    
                   double* bCb, double* d_invD_d, double* log_det_Chol, 
		   double* log_det_invD, WORK* myWork);



int accept_new_lambda(QTL_INFO* aqtl, QTL_INFO** all_qtls, int nn, 
					  int newMark, double newPos, int gmiss,
                      params* theparams, int bc, int revjump);
void Gibbs_update_geno(QTL_INFO* aqtl, int nn, 
			 	  	   double** p_resid, double sigmasq,  
                       params* theparams, int bc, int gmiss,
					   double* new_r, 
					   int** p_genotype, double** p_newResid,
					   double* log_norm_const, int revjump);
void propose_qtl_geno(QTL_INFO* aqtl, int nn, 
			 		 double* resid, double sigmasq,  
                     params* theparams, int bc, int gmiss,
					 double* new_r, int* genotype, double* newResid,
					 double* log_norm_const, int revjump);
void genProbs(int nn, params* theparams, int bc, int gmiss, 
			  QTL_INFO* aqtl, int revjump);

void fixed_locus_update(int nn, int nqtl, QTL_INFO** all_qtls,
			MCMC_PARAM* myMCMC, DATA* myData, PRIORS* priors, 
			WORK* myWork,
			params* theparams, double** nmove);

int birth(DATA* myData, CHROMOSOME* chrom, PRIORS* priors, 
	  params* theparams, WORK* myWork,
	  MCMC_PARAM* myMCMC);

void initQtl(int nn, QTL_INFO* aqtl, 
	     params* theparams, int bc, int gmiss,
	     int revjump);
void setupChromosomes(DATA* myData, MCMC_PARAM* myMCMC, CHROMOSOME** pChromInfo,
		      params* theparams, int* markers, 
			  int *chrarray, int *nmararray, double *posarray, double *distarray,
			  int* markersPerChrom);
void setupMCMC(MCMC_PARAM* myMCMC, DATA* myData, WORK* myWork, PRIORS* priors,
	       params* theparams, int* genodata,
    	   int *chrarray, int *nmararray, double *posarray, double *distarray,
		   int* markersPerChrom, double* pheno, CHROMOSOME** pChromInfo,
		   int* revjump,int* nind,double* burn_in,
	       double* preburn_in,int* niter,int* nby, double* mu, double* sigmasq,
	       double* meanmu, double* sdmu, double* siga1, double* siga2,
	       double* meanadd, double *sdadd,double* meandom, double*sddom,
	       double* qtlmean);
void setupTraitData(DATA* myData, int nn, double* traity);
void setupQtl(DATA* myData, MCMC_PARAM* myMCMC, WORK* myWork, 
	      CHROMOSOME* chromInfo, PRIORS* priors, 
	      int* new_chrom, double* new_lambda, 
	      params* theparams);			  
#endif

