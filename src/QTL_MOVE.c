/**************************************************************
    File:       main.c
    Written by: Patrick Gaffney
    Date:       November 11, 2000
    Version:    0.4

    Purpose:

      To move the QTL genome records in the genome linked
	list supplied in QTLCart.

************************************************************/

#include "revjump.h"
#include "ranlib.h"

/****************************************************************************
 These routines utilize the genome structure provided in QtlCart
 (which stores position of markers and inter-marker distances).
 QTLs genome records are stored permanently in the genome structure
 (position and intermarker distance).  So when calc_cond_prob is 
 called with the qtl genome record, it returns the probability of 
 possible QTL genotype (AA,Aa and aa) conditional on the available 
 flanking marker (or QTL) data).

 The two routines used are:
  (1) insertQtl
  (2) removeQtl

 Both require a pointer to the qtl of interest.  The former also 
 requires a reference pointer to ANY elelemt in the genome list 
 (ideally close to the insertion point), the lambda position 
 (xx.yyy where xx=leftmost marker # and .yyy is the fraction of
 the distance to the next marker).  Also required are the chromosome
 number and the qtl number. 

 ***************************************************************************/



double calc_intermarker_dist(CHROMOSOME* chrom, int m)
/* calc the distance from marker m to marker m + 1 */
{
   mygenome* lmark, *rmark;
   
   assert(m >= 1 && m <= chrom->nMark);

   lmark = chrom->mark_genome[m];
   rmark = chrom->mark_genome[m+1];
   return rmark->pos - lmark->pos;
}


int EQUALS(double x, double y)
{
   return (ABS(x-y) < EPSILON);
}


int LE(double x, double y)
{
   return (x-y <= EPSILON);
}

int LT(double x, double y)
{
   return (x-y < -EPSILON);
}


double getPos(CHROMOSOME* chrom, double lambda, mygenome** lptr)
{
  /* converts from Jaya's coding for lambda [xx.yyyy where xx gives 
     the number of the nearest left marker (1..nmark) and .yyyy
     (0 to 1.0-tol) is the fraction of the inter-marker distance 
     from said left marker] to an absolute position (in Morgans)
     along the chromosome.  

     Output:
        converted value
        lptr pointes to leftmost marker


  */
  mygenome* rptr;

  *lptr = chrom->mark_genome[(int)lambda];
  rptr = chrom->mark_genome[(int)lambda+1];

  return (*lptr)->pos + (rptr->pos - (*lptr)->pos)*(lambda - (int)lambda);
}




void insertQtl(mygenome* qptr, int lmark, double pos, CHROMOSOME* chrom, int qtlNum)
{
  /*  inserts Qtl genome record.  qptr can be blank.  All details will be added.
	  qtlNum MUST be positive (>0), and should be the index of the (ordered) QTL along chromosome 
	*/
        
     mygenome* genome_lmark;
     mygenome* lptr;   /* may be the same as lmark */
     double old_dist;   
 
     qptr->chrom = chrom->num;
     qptr->markr = -qtlNum;
	 genome_lmark = chrom->mark_genome[lmark];

     qptr->pos = pos;
     lptr = findInsertPos(qptr->pos, genome_lmark, qptr->chrom);

	qptr->prev = lptr;
	qptr->next = lptr? lptr->next: NULL;

	if (qptr->next) qptr->next->prev = qptr;
	if (lptr) lptr->next = qptr;   /* recall, qptr->prev is lptr */

	/* next the dist/pos */
        old_dist = lptr->dist;
	lptr->dist = qptr->pos - lptr->pos;
	qptr->dist = old_dist - lptr->dist;
}


void removeQtl(mygenome* qptr)
{
	mygenome* prev, *next;

	prev = qptr->prev;
	next = qptr->next;

   /* first adjust the prev */
   if (prev) {
      prev->next = next;
      prev->dist = (next)? next->pos - prev->pos: 0.0;
	              /* change if tail DNA */
   }

   /* now the next marker/qtl */
   if (next) next->prev = prev;	

}


int checkIntegrity(int nQtl, CHROMOSOME* chrom)
{
   mygenome* g;
   int chrNum;
   int mark, qtl;

   if (chrom->nMark == 0) return 1;

   g = chrom->mark_genome[1];  /* first marker */
   chrNum = g->chrom;
   
   if (!((g->pos == 0.0 && g->markr == 1) || g->markr <=0 && chrNum == chrom->num))
	   printf("hello");

   assert((g->pos == 0.0 && g->markr == 1) || g->markr <=0 && chrNum == chrom->num); 
   g = g->next;

   /* check that distances are equal */
   for (mark=1,qtl=0; g && mark+qtl < chrom->nMark + chrom->nQtl; 
                      g = g->next)
   {
     if (!EQUALS(g->pos - g->prev->pos, g->prev->dist)) 
        printf("Error in distances %d %d %g\n", mark, qtl, g->pos); 
	 else if (g->chrom != chrNum)
        printf("Error in chrNum, markr %d chr %d != %d\n", mark, qtl, g->pos); 
	 else if (g->markr > chrom->nMark || g->markr < -nQtl)
        printf("Error, marker num (%d) should be %d ..  %d \n", g->markr, -nQtl, chrom->nMark); 

     if (g->markr >=0) mark++;
     else qtl++;    
   }

   if (mark != chrom->nMark || qtl != chrom->nQtl)
	   printf("Error, number mark %d != chrom nMark %d, or nQtl %d != chrom nQtl %d\n",
	           mark,chrom->nMark,qtl,chrom->nQtl);
   return 1;
}


void replaceQtl(mygenome* old, mygenome* new)
{
   new->dist = old->dist;
   new->pos = old->pos;
   new->prev = old->prev;
   new->next = old->next;
   restoreQtl(new);
} 



/* assumes pos, prev and next fields are correct */
void restoreQtl(mygenome* qptr)
{
	mygenome* prev, *next;

	prev = qptr->prev;
	next = qptr->next;

  if (prev) {
    prev->next = qptr;
    prev->dist = qptr->pos - prev->pos;
	if (prev->dist <0) 
		printf("<* (%d) %g (%d) %g **\n", qptr->markr, qptr->pos, qptr->prev->markr, qptr->prev->pos);
	
    assert(prev->dist >=0);
  }
  
  if (next) {
    next->prev = qptr;
    qptr->dist = next->pos - qptr->pos;

	if (qptr->dist <0) 
		printf("++ (%d) %g (%d) %g +>\n", qptr->next->markr, qptr->next->pos, qptr->markr, qptr->pos);
    
	assert(qptr->dist >=0);
  }
  else
    qptr->dist = 0;
}




mygenome* findInsertPos(double pos, mygenome* lmark, int chrom)
{
  /*      We ASSUME that there is always a marker to the left of the qtl
     to be inserted (lmark).  If ever you intend to insert Qtls in the tail
     DNA, outside known markers, it is necessary to add
     non-real markers with markr set to 0.  The individuals' records
     must be updated so that these 'markers' have a missing genotype
     value (non-trivial amount of work).
  */
   mygenome* next;

   while ((next=lmark->next) && next->pos < pos && next->chrom == chrom) 
	   if (next->next == lmark->prev) 
	   {
		   break;
	   }
	   else 
         lmark = next;
   
   /* one final check */
    if (!(lmark && chrom == lmark->chrom))
	   printf("hello2");
   assert (lmark && chrom == lmark->chrom);
   return lmark;
}







