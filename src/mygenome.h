#ifndef MYGENOME_GUARD
#define MYGENOME_GUARD

typedef struct MyGenome { /*Structure to hold a genome defined by the genetic linkage map*/
  int chrom;                 /* Chromosome of the marker interval */
  int markr;                 /* The marker that precedes the interval, can be 0 if borders are allowed */
  double dist;               /* This distance, from marker markr to markr+1 on chromosome chrm, is in M */
  double pos;                /* Position of marker from left telomere in Morgans*/
  int* genotype;             /* marker genotype */
  struct MyGenome *prev;            /* Pointer to previous marker interval */
  struct MyGenome *next;             /* Pointer to next marker interval */
}  mygenome;



#endif

