/**************************************************
  File:         update_parms.c
  Written by:   Patrick Gaffney
  Date:         November 11, 2000
  Version:      0.4

  Purpose:

    Update the model parameters.
**************************************************/

#include "revjump.h"
#include "ranlib.h"
#include "qtl.h"

/* ------------------------------------------------------------------
   Functions in this file: 

   update_effects:   proposes updated effects

* ------------------------------------------------------------------*/

void update_effects(int nn, int nQtl, int revjump, double* y, 
		    double* mu, double sigmasq, int* na, 
		    QTL_INFO** qtls, PRIORS* priors,WORK* myWork)
     /* 
   A multivariate Gibbs Sampler for mean, additive and dominance effects. Uses
   Cholesky decomposition.  Underlying method described in generate_effects().
 

  Inputs
  ======
  nn         the number of plants
  nQtl       number of QTL in the model
  y          vector 1,..,nn of phenotypes
  sigmasq    environmental variance
  na         number of effects of each type in the model with nQtl QTL
  w          vector of weights, indexed by 1..na[ADD]+na[DOM]+1. The first weight should be 1. 
             The remainder correspond to weights for additive/dominance effect variances of 
			 each QTL.   Thus if the first QTL had an additive and dominance effect, w[2] 
			 would give the additive effects variance weight and w[3] that for the domiance 
			 effect,...
  qtls       information for the nQtl QTL (1..nQtl)
  priors     prior structure (means and variance for effects and mean)
  myWork     work structure (be aware of side-effects)

  Outputs
  =======
  mod_effect the result (generally part of myWork, so again beware).  Indexed 1,..,na[ADD]+
             na[DOM]+1.  The first term is the proposed mean.  Additive and dominance effects
			 for each QTL in turn (all effects for one QTL followed by those for the next QTL)
			 are stored.


  Side-effects
  ============
      myWork ->  u, k, v. wr
	  
       
*/
{
  /* set up work variables */
  double* pmean = myWork->pmean;
  double* pvar = myWork->pvar;
  double** chol = myWork->chol;
  double* p = myWork->p;
  double** XtX = myWork->XtX;
  double* XtY = myWork->XtY;
  double* mod_effect = myWork->mod_effect;
  int nterm;
  
  double bCb, d_invD_d;
  double log_det_Chol, log_det_invD;

  setPriorMeanVar(nQtl, revjump, qtls, NULL,
		  priors, myWork->weight, na, pmean,pvar);
  nterm = setAddDomCovMatrix(nn, nQtl, qtls, y, XtX, XtY, na);

  if (na[ADD] + na[DOM] + 1 != nterm) {
    printf("error in na in update_effects %d + %d != %d\n",na[ADD], na[DOM], nterm);	         
    exit(1001);
  }

  generate_effects(nterm, XtX, XtY, pmean, pvar, 
                   sigmasq, chol, p, myWork->u, mod_effect,
		   &bCb, &d_invD_d, &log_det_Chol, &log_det_invD, 1);

  setEffect(nn, nQtl, y, myWork->mod_effect, qtls, 
	    mu, myWork->weight, myWork->resid, na);
  setCholParams(myWork, bCb, d_invD_d, log_det_Chol, log_det_invD);
}







void setPriorMeanVar(int nQtl, int revjump, QTL_INFO** qtls, QTL_INFO* newQtl,
		     PRIORS* priors, 
		     double* w, int* na, double* pmean, double* pvar)
     /* setup prior means & variances, and the weights for the model with nQtl QTL */
{
  int i, idx;
  QTL_INFO* aqtl;

  pmean[1] = priors->mean[MU];
  pvar[1] = priors->var[MU];  
  w[1]=1;
  setWeights(nQtl, na, revjump, qtls, w, newQtl);

  for (i=1, idx = 2; i<=nQtl-1; i++)
    {
      aqtl=qtls[i];
      if (aqtl->flag & QTL_ADD) {
	pmean[idx] = priors->mean[ADD]; 
	pvar[idx] = priors->var[ADD]*w[idx]; 
	idx++;
      }
      if (aqtl->flag & QTL_DOM) 
	{
	  pmean[idx] = priors->mean[DOM]; 
	  pvar[idx] = priors->var[DOM]*w[idx]; 
	  idx++;
	}
    }

  /* now take care of the newQtl */
  if (!newQtl && nQtl > 0) newQtl = qtls[nQtl];
  if (newQtl)
    {
      if (newQtl->flag & QTL_ADD) {	
	pmean[idx] = priors->mean[ADD]; 
	pvar[idx] = priors->var[ADD]*w[idx]; 
	idx++;
      }
      if (newQtl->flag & QTL_DOM) 
	{
	  pmean[idx] = priors->mean[DOM]; 
	  pvar[idx] = priors->var[DOM]*w[idx]; 
	  idx++;
	}
    }

}





int generate_effects(int nterm, double** XtX, double* XtY, 
		     double* pmean, double* pvar,double sigma2, 
		     double** chol, double* p, double* u, double* mod_effect,
		     double* bCb, double* d_invD_d,
		     double* log_det_Chol, double* log_det_invD, int genBeta)

{
  /* 
     Generates values from a multivariate normal with mean c and covariance inv(C).
     We don't specify C and c directly to this routine.  The following paragraphs
     describe how they are computed.  

  
     Input:
     start        set to 1 to compute cholesky, 0 or less to prevent computation of cholesky,
     and to a value less than or equal to nterm to recompute these elements of
     the cholesky.
     nterm        number of terms in the s2A matrix
     s2A
     s2Aa
     k
     v
     s2
     myWork
     
     Output:
     myWork->add
     myWork->u


   
	 
     In Bayes analysis, we oftentimes have a likelihood which has a normal density, 
     and a prior (congugate) which is also normal.  This leads to a normal with 
     exponent of the following form:

     (u-a)'A(u-a) + (u-d)'inv(D)(u-d) 
					  
     This equals (less a constant wtr u) 
            
     (u-c)'C(u-c) 

     where c=inv(A+B)(Aa+inv(D)d) and C = A+inv(D).  In our case D is diagonal, 
     with D = diag(pvar), and d=pmean.  pmean and pvar are the prior additive 
     effect mean and variances respectively.

     We are given (as inputs to this fuction): XtX = (s^2)A and XtY = (s^2)Aa, 
     pmean = d and pvar = diag(D), where s^2 is the environmental variance.
   
     Let M = (s^2)C.

     Then M = (s^2)A + (s^2)inv(D) = XtX + (s^2)inv(D) and 
     c=inv(M)(XtY + (s^2)inv(D)d) = inv(M) x,  or equivalently if M is 
     of full rank (which it is), then Mc = x and we solve by obtaining Cholesky 
     for M and employing a forward-backward substitution to solve for c.  
   
     To generate the multi-variate normal, we need z = s * choleskyL(inv(M)) * u, 
     where u is multivariate normal with mean 0 and variance I and choleskyL() 
     returns the lower triangular matrix from the Cholesky decomposition, and s is
     the standard deviation (sqrt of environ variance). Since choleskyL(M) * z = s*u
     we can solve for z using forward substitution.


 */
  double t;
  int i, info;
  double s = sqrt(sigma2);

  /* see Chapter 4 of Thesis, proposing birth effects (not Jaya's) for
     an explanation of these terms. */
  *bCb = 0.0;
  *d_invD_d = 0.0;
  *log_det_Chol=0.0;  
  *log_det_invD = 0.0;


  /* effectively make s2A a p.d. matrix, we merge likelihood and prior [both
     normal] to form the posterior for effects  */
  for (i=1; i<=nterm;i++) 
    {
      t = sigma2/pvar[i]; 
      XtX[i][i] += t;  
      XtY[i] += pmean[i]*t; 
      *log_det_invD -= log(pvar[i]);
      *d_invD_d += (pmean[i] * pmean[i])/pvar[i];
    }

  /* now we form Cholesky if needed */
  mycopy(nterm, nterm, chol, XtX);                  /* we don't want to destroy XtX  */
  info = cholesky(nterm, chol, p);                  /* find L s.t. XtX = LL'         */
	                                
  /* find mean of MVN for effects posterior */
  choleskySolve(nterm,chol,p, XtY, mod_effect);     
                                                          

  /* needed for acceptance term in birth/death */
  for (i=1; i<=nterm; i++) *log_det_Chol += log(p[i]);
  for (i=1; i<=nterm; i++) *bCb += XtY[i] * mod_effect[i];
  *bCb /= sigma2;

  if (genBeta)                /* if we want to simply regenerate Cholesky, genBeta==0 */                
    {  
      for (i=1; i<=nterm; i++) u[i] = s * snorm();     /* generate random numbers       */
      backwardSubs(nterm, chol, p, u);                  /* forward substitution to get z */                                                           
      for (i=1; i<=nterm; i++) mod_effect[i] += u[i];  /* add to mean of MVN            */
    }

  /* restore the XtX matrix and the XtY vector to their original forms */
  for (i=1; i<=nterm;i++) {t = sigma2/pvar[i]; XtX[i][i] -= t;  XtY[i] -= pmean[i]*t;}

  if (info != 0)
    printXtX(nterm,XtX);

  return info;
}





  



void choleskySolve(int n, double **a, double* p, double* b, double* x)
     /* returns x where Ax=b  and a,p are the cholesky decomposition of A   */
     /* neither a nor b is destroyed                                        */
{
  int i,k;
  double sum;

  for (i=1;i<=n;i++) {
    for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
    x[i]=sum/p[i];
  }
  for (i=n;i>=1;i--) {
    for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }
}


void forwardSubs(int n, double **a, double*p, double* b)
     /* assumes a is lower triangular, with its diagonal in p, b is destroyed */
{
  int i,j;
  double sum;

  for (i=1;i<=n;i++) 
    {
      sum=b[i];
      for (j=1;j<=i-1;j++) sum -= a[i][j]*b[j];
      b[i]=sum/p[i];
    }
}


void backwardSubs(int n, double **a, double*p, double* b)
     /* assumes a is upper triangular, with its diagonal in p, b is destroyed */
{
  int i,k;
  double sum;

  for (i=n;i>=1;i--) 
    {
      sum=b[i];
      for (k=i+1;k<=n;k++) sum -= a[k][i]*b[k];
      b[i]=sum/p[i];
    }
}


int cholesky(int n, double **a, double* p)
     /* this destroys the original a matrix. since now it contains the cholesky decomposition 
   If the matrix is singular, we have problems */
{
  int i,j,k;
  double sum;

  for (i=1;i<=n;i++) {
    for (j=i;j<=n;j++) {
      for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
      if (i == j) {
	if (sum <= 0.0)
	  return i;
	p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }
  return 0;
}


double log_determinant(int nterm, double* p, double sigmasq)
{
  int i;
  double det=0;
  double log_d;
  static double logTwoPi = 1.837877066;

  log_d = 0.5 * (log(sigmasq) + logTwoPi);

  for (i=1; i<=nterm; i++) det += log(p[i]);

  det -= nterm*log_d;
  return det;
}

int incremental_cholesky(int n,  double** a, double* p, double * newcol)
     /* chol and p contain the results of a Cholesky decomposition of an n by n matrix.  We
   now wish to obtain the decomposition of the original symmetric matrix with a new
   column/row added (newcol) */
{
  int i,k;
  double sum;

  for (i=1;i<=n;i++) 
    {
      sum = a[i][n+1] = newcol[i];
      for (k=i-1;k>=1;k--) sum -= a[i][k]*a[n+1][k];
      a[n+1][i] = sum/p[i];
    }

  for (sum=newcol[n+1],k=n;k>=1;k--) sum -= a[n+1][k]*a[n+1][k];
  if (sum <=0) return n+1;
  else p[n+1] = sqrt(sum);

  return 0;
}




void mycopy(int n, int m, double** to, double** from)
{
  int i, j;

  for (i=1; i<=n; i++)
    for (j=1; j<=m; j++)
      to[i][j] = from[i][j];
}






int setAddDomCovMatrix(int nn, int nQtl, QTL_INFO** selectedQtls, double* y,
		       double** XtX, double* XtY, int* na)
{
  /*  Evaluates matrix add[i][j] = sum(x[i][k] * x[j][k], k=1,..,nn) where x[i][k] 
      is the genotype of the ith QTL for the kth individual.  
	
      Also calculates sum(x[i][k] * r[k], k=1,..nn) where  
	 
      r[k]=resid[k]           if residIsY != 0 
      and
      r[k] = resid[k] - sum(selectedQtls[j]->a[ADD] * x[i][k], j=1,..,nQtl)  otherwise
		 
      i.e. r[k] is the residual when the selectedQtls are not included in the model.

      We assume selectedQtls have the flag & QTL_ADD evaluating to true (i.e. the QTL
      has an additive effect 
  */

  int i;
  QTL_INFO* aqtl;
  int old_na[3];
  int type;
	

  /* calculate the row 1 (mean-effect) terms and sum of genotype * phenotype */			
  XtX[1][1]=nn;
  for (i=1, XtY[1]=0.0; i<=nn; i++) XtY[1] += y[i];
  old_na[ADD]=old_na[DOM]=0;

  /* now calculate off-diagonal terms */
  for (i=1; i<=nQtl; i++)
    {
      aqtl = selectedQtls[i];
      addColToAddDom(nn, i-1, old_na, selectedQtls, aqtl,
		     y, XtY, XtX);

      /* we need to adjust old_na for addColToAddDom */
      for (type=ADD; type<=DOM; type++) 
	if (aqtl->flag & type) old_na[type]++;
    }

  return old_na[ADD]+old_na[DOM]+1;   /* dimension of matrix */
}




void setAddDomDiag_Row1(int nn, int idx, QTL_INFO* aqtl, double* y, double** XtX, double* XtY)
{
  int k;
  double x,z;
  int* geno;
  double sum_x, sumSqr_x, wr_x;
  double sum_z, sumSqr_z, wr_z, wr_zx;

  geno = igenotype(aqtl);

  if (aqtl->flag & QTL_ADD)
    {
      sum_x= sumSqr_x = wr_x = 0.0;

      for (k=1; k<=nn; k++) {
	x = geno[k]; 
	sum_x += x; sumSqr_x += (x*x); wr_x += (y[k]*x);
      }

      if (XtX) 
	{
	  XtX[1][idx] = XtX[idx][1] = sum_x;
	  XtX[idx][idx] = sumSqr_x;
	}
      XtY[idx] = wr_x;
      idx++;
    }

  if (aqtl->flag & QTL_DOM)
    {
      sum_z= sumSqr_z = wr_z = wr_zx = 0.0;

      for (k=1; k<=nn; k++) {
	z = DOM_GENO_VALUE(geno[k]); 
	sum_z += z; sumSqr_z += (z*z); wr_z += (y[k]*z); wr_zx += z*geno[k];
      }

      if (XtX) 
	{
	  XtX[1][idx] = XtX[idx][1] = sum_z;
	  XtX[idx][idx] = sumSqr_z;
	  if (aqtl->flag & QTL_ADD) 
	    XtX[idx][idx-1] = XtX[idx-1][idx] = wr_zx;
	}
      XtY[idx] = wr_z;
    }
}


void addColToAddDom(int nn, int old_nQtl, int* old_na, QTL_INFO** qtls, QTL_INFO* aqtl,
		    double* y, double* XtY, double** XtX)

     /*
   Adds one or more columns to the XtX matrix for the new QTL.  Currently only columns
   for the additive and dominance effects are considered.  Also adds the corresponding
   entry(s) in the XtY matrix.


   Inputs
   ======
    nn      ... number of individuals
	old_nQtl .. number of QTL already in the model
    old_na  ... number of additive and dominance effects in model before we wish to add
	            aqtl
	qtls    ... array of pointers to QTL record
	aqtl    ... pointer to QTL record to be added
	y       ... phenotype vector
	XtY     ... transpose of design matrix times the phenotype vector (valid for 
	            old_nQtl model)
	XtX     ... for old_nQtl model


   Outputs
   =======
    XtY     ... now valid for model with old_nQtl+1 QTL
	XtX     ... now valid for model with old_nQtl+1 QTL

*/    
{
  int i,idx,idx2;
  QTL_INFO* bqtl;

  idx = 2 + old_na[ADD] + old_na[DOM];
  setAddDomDiag_Row1(nn, idx, aqtl, y, XtX, XtY);

  for (i=1, idx2=2; i<=old_nQtl; i++)
    {
      bqtl = qtls[i];
      calc_interaction(idx, idx2, XtX, aqtl, bqtl, nn);
      idx2 += bqtl->nParam;
    }
}


void calc_interaction(int idx, int idx2, double** a, QTL_INFO* aqtl, QTL_INFO* bqtl, int nn)
{
  /* aqtl is associated with index idx, bqtl with index idx2 */
  double sum;
  int type, type2;
  int k;
  int* geno1, *geno2;
  int copy_idx2=idx2;


  geno1 = igenotype(aqtl);
  geno2 = igenotype(bqtl);

  for (type=ADD; type<=DOM; type++)
    if (aqtl->flag & type)
      {
	for (type2=ADD, idx2 =copy_idx2; type2<=DOM; type2++)
	  if (bqtl->flag & type2)
	    {
	      if (type ==ADD && type2 == ADD)
		{
		  for (k=1, sum=0; k<=nn; k++) sum += geno1[k]*geno2[k];
		}
	      else if (type==ADD)
		{
		  for (k=1, sum=0; k<=nn; k++) sum += geno1[k]*DOM_GENO_VALUE(geno2[k]);
		}
	      else if (type2==ADD)
		for (k=1, sum=0; k<=nn; k++) sum += geno2[k]*DOM_GENO_VALUE(geno1[k]);
	      else
		for (k=1, sum=0; k<=nn; k++) sum += DOM_GENO_VALUE(geno1[k])*DOM_GENO_VALUE(geno2[k]); 	              
	      a[idx][idx2] = a[idx2][idx] = sum;
	      idx2++;
	    }
	idx++;
      }
}




void removeRowCol(int n, double** a, int row)
{
  int i,j;

  for (i=1; i<row; i++) 
    for (j=row+1; j<=n; j++)
      a[i][j-1] = a[i][j];

  for (i=row+1; i<=n; i++) if (i!=row)
    SwapDVec(&a[i], &a[i-1]);
}


void moveQtlToEndOfXtX(int nQtl, QTL_INFO** qtls, int changeQtl, double** XtX, double* XtY, 
		       int nterm)
{
  int i,j,k;
  int idx, idx2;
  QTL_INFO* cqtl = qtls[changeQtl];
  double* z[2];
  double zz[2];
  int nParam = cqtl->nParam;

  for (i=1, idx=2; i<changeQtl; i++) idx += qtls[i]->nParam;

  /* swap rows, move remaining rows to fill void left by moving QTL */
  for (idx2=idx, j=0; j< nParam; j++, idx2++) z[j] = XtX[idx2];
  for (j=idx, k=idx2; k<=nterm; j++, k++) XtX[j] = XtX[k];
  for (j=0; j< nParam; j++) XtX[nterm-j] = z[nParam-j-1];

  /* now clean up the columns */
  for (i=1; i<=nterm; i++)
    {
      for (j=0; j< nParam; j++) zz[j] = XtX[i][j+idx];
      for (j=idx, k=idx2; k<=nterm; j++, k++) XtX[i][j] = XtX[i][k];
      for (j=0; j<nParam; j++) XtX[i][nterm-j] = zz[nParam-j-1];
    }

  /* now exchange the XtY */
  for (j=0; j< nParam; j++) zz[j] = XtY[j+idx];
  for (j=idx, k=idx2; k<=nterm; j++, k++) XtY[j] = XtY[k];
  for (j=0; j< nParam; j++) XtY[nterm-j] = zz[nParam-j-1];

  for (i=changeQtl+1; i<=nQtl; i++) {
    qtls[i-1] = qtls[i];
    qtls[i-1]->qptr->markr = -(i-1);
  }
  qtls[nQtl] = cqtl;
  cqtl->qptr->markr = -nQtl;
}




void swapRowAndCol(int nn, double** a, int row1, int row2)
{
  int i;

  /* interchanges two row-column combinations */
  for (i=1; i<=nn; i++)
    SwapDble(&a[row1][i], &a[row2][i]);

  for (i=1; i<=nn; i++)
    SwapDble(&a[i][row1], &a[i][row2]);
}






/*====================================================================

  Gibbs Sampler for mean, additive and dominance effects for one
  parameter at a time

 ====================================================================*/


double update_mu(int nn, double mu, double sigmasq, double* resid,
		 double mu_prior_mn, double mu_prior_var)
     /* ------------------------------------------------------------------
    simple Gibbs sampler for mean
 * ------------------------------------------------------------------*/
{
  double mean0, var0, mu_val;
  int i;
  double total, diff;

  total = 0.0;
  for(i=1; i<=nn; i++) {
    total += (resid[i] + mu);
  }

  var0 = (1.0/mu_prior_var) + nn/sigmasq;

  mean0 = (mu_prior_mn/mu_prior_var + total/sigmasq) / var0;
  mu_val = mean0 + gennor(0.0,1.0)/ sqrt(var0);

  diff = mu - mu_val;
  for(i=1; i<=nn; i++) resid[i] += diff;  
  if (mu_val < -100 || mu_val > 100)
    return mu_val;
  return mu_val;
}




double update_add_effect(int nn, int nQtl, double sigmasq,
			 QTL_INFO* modQtl, double* resid, 
			 double effect_prior_mn, double effect_prior_var)
     /* ------------------------------------------------------------------
    simple Gibbs sampler for additive effect  for QTL given by modQtl
	Updates the residual.
 * ------------------------------------------------------------------*/
{
  double work1, work2, work3, var1, mean1, new_a, old_a, u;
  int i;
  int* geno = igenotype(modQtl); 

  old_a = modQtl->a[QTL_ADD];
	
  work1=work3=0.0;
  for(i=1; i<= nn;i++){ 
    if (*geno)
      {
	work2 = resid[i] + old_a * geno[i];
	work1 += (geno[i] * work2);
	work3 += (geno[i] * geno[i]);
      }
  }

  var1 = 1/effect_prior_var + (work3/sigmasq);
  mean1 = (effect_prior_mn/effect_prior_var + work1/sigmasq) / var1;
    
  u = gennor(0.0, 1.0);
  new_a = mean1 + u / sqrt(var1);

  /* change residuals */
  for (i=1; i<= nn; i++) 
    resid[i] -= (new_a - old_a) * geno[i];

  return new_a;
}



double update_dom_effect(int nn, int nQtl, double sigmasq,
			 QTL_INFO* modQtl, double* resid, 
			 double effect_prior_mn, double effect_prior_var)
     /* ------------------------------------------------------------------
    simple Gibbs sampler for dominance effect for QTL given by modQtl
	Updates the residual.
 * ------------------------------------------------------------------*/
{
  double work1, work3, var1, mean1;
  int i;
  int* geno = igenotype(modQtl);
  double old_d, new_d;

  old_d = modQtl->a[QTL_DOM];

  work1=work3=0.0;
  for(i=1; i<= nn;i++){ 
    if (geno[i]==0)
      {
	work1 += (resid[i] + old_d);
	work3 += 1; 
      }
  }

  var1 = 1/effect_prior_var + (work3/sigmasq);
  mean1 = (effect_prior_mn/effect_prior_var + work1/sigmasq) / var1;
    
  new_d = mean1 + gennor(0.0, 1.0) / sqrt(var1);

  /* change residuals */
  for (i=1; i<= nn; i++) 
    if (geno[i]==0)
      resid[i] -= (new_d - old_d);

  return new_d;
}


/* ------------------------------------------------------------------
    old simple Gibbs sampler for variance
 * ------------------------------------------------------------------*/

double Gibbs_update_sigmasq(int nn, double* resid,
			    double prior_sigsq_a1,
			    double prior_sigsq_a2)
{
  double sum, work, sigmasq_val;
  int i;

  sigmasq_val = gengam(1.0, prior_sigsq_a1);
  sum = 0.0;

  for (i=1; i<=nn; i++) 
    {
      work = resid[i];
      work *= work;
      sum += work;
    }

  sigmasq_val = (prior_sigsq_a2 + sum/2.0) / sigmasq_val;

  return sigmasq_val;
}


double MH_update_sigmasq(int nn, double* resid,
			 double sigmasq,
			 double deltaS2, double maxS2)
{  
  double residSS, work, new_sigmasq;
  double log_ratio, uni_ran;
  int i;

  /* calculate residual SS */
  for (i=1, residSS=0.0; i<=nn; i++) 
    {
      work = resid[i];
      work *= work;
      residSS += work;
    }
  new_sigmasq = sigmasq + genunf(-1,1) * deltaS2;
  if (new_sigmasq > maxS2) new_sigmasq = 2*new_sigmasq - maxS2;
  else if (new_sigmasq < 0) new_sigmasq = -new_sigmasq; 
  
  log_ratio = nn/2.0 * log(sigmasq/new_sigmasq) - 0.5 * residSS * (1/new_sigmasq - 1/sigmasq);
  uni_ran = log(genunf(0.0,1.0));
  if (uni_ran < log_ratio) return new_sigmasq;      /* accept proposal */
  else return sigmasq;
}

double calc_prior_effect_term(double effect, double w, double effect_mean,
			      double log_effect_var, double effect_var)
{
  static double logTwoPi = 1.837877066;
  double temp;

  temp = effect - effect_mean;
  temp *= temp;
  temp /= -(2 * w * effect_var);
  temp -= 0.5*(logTwoPi + log(w) + log_effect_var);

  return temp;
}



double qtl_prior_effect_ratio (double effect_prop, double w_prop,
			       double effect, double w,
			       double effect_mean, double effect_var, 
			       double log_effect_var)
{
  double temp;

  if (w_prop == w)
    {
      temp = effect_prop - effect;
      temp *= (effect_prop + effect - 2*effect_mean);
      temp /= -(2 * effect_var * w);
    }
  else if (effect == effect_prop)
    {
      temp = (effect_prop - effect_mean);
      temp = -(1/w_prop - 1/w) * temp * temp / (2 * effect_var);
      temp -= 0.5 * log(w_prop/w);
    }
  else
    {
      temp = calc_prior_effect_term(effect_prop, w_prop, effect_mean, 0.0, effect_var) -
	calc_prior_effect_term(effect, w, effect_mean, 0.0, effect_var);
    }

  return temp;
}




double update_effect_prior_var(int type, int nQtl,  QTL_INFO** all_qtls,
			       PRIORS* priors, double y_var)
     /* ------------------------------------------------------------------
   Updates the prior variance for either ADD or DOM effects (specifed
   by type) for a model with nQtl QTL.
 * ------------------------------------------------------------------*/
{
  double u, old_var, new_var;
  double effect, weight;
  double old_r, new_r, ratio_tau, ratio_accept;
  double DELTA = 0.1;
  int j;
  double t;
  QTL_INFO* aqtl;
	
  if (!priors->sampleVar[type]) return -1;

  old_var = priors->var[type];
	                           
  t = old_var/2.0/y_var;
  old_r = pow(t, priors->alpha[type]-1) * 
    pow(1.0 - t, priors->beta[type]-1);
  u = genunf(-1,1) * DELTA * y_var;
  new_var = old_var + u;

  if (new_var < 0) new_var = -new_var;
  if (new_var == 0) new_var = EPSILON;
  else if (new_var > 2*y_var) new_var = 4*y_var - new_var;
	

  t = new_var/2.0/y_var;
  new_r = pow(t, priors->alpha[type]-1) *
    pow(1.0 - t, priors->beta[type]-1);


  ratio_tau = new_var/old_var;
  ratio_accept = 0;

  for (j=1; j<=nQtl; j++)
    {
      aqtl = all_qtls[j];
      if (type & aqtl->flag) 
	{
	  effect = aqtl->a[type];
	  weight = aqtl->w[type];
	  ratio_accept += 
	    qtl_prior_effect_ratio ( effect, weight * ratio_tau, effect, weight,
				     priors->mean[type], old_var,priors->var[type]);		
	}
    }
  ratio_accept = exp(ratio_accept) * (new_r/old_r);


  u = RANF();

  if (u < ratio_accept)  /* make move */
    {
      priors->var[type] = new_var;
      priors->log_var[type] = log(new_var);
    }		

  return priors->var[type];
}


