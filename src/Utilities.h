#ifndef UTILITIES_GUARD
#define UTILITIES_GUARD
/*  
          Utilities.h
          
Copyright (C) 1996 Christopher J. Basten, Bruce S. Weir and Zhao-Bang Zeng.

This file is part of QTL Cartographer. QTL Cartographer is free software; you
can redistribute it and/or modify it under the terms of the GNU  General
Public License as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

QTL Cartographer is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along with
QTL Cartographer; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

extern int whosemf;
extern long     ix, Ix;
extern double mapparam;


extern FILE *fileopen(char *name, char *mode);  /* open a file with error messages */
extern void fileclose(char *name, FILE *fp);  /* close a file with error messages */

extern void MemoryAccount(char *buffer);
extern double lrtolod(double lr);
extern double lodtolr(double lod);
extern int is_number(char *buffer);
extern int is_pinteger(char *buffer);
extern int myfgets(char *buffer, int nn, FILE *fptr);

extern void shuffle_ivector(int *v, int lb, int ub);
extern double mapfunc(double value, int flag);
extern double iKosambi(double mm); 
extern double Kosambi(double rr); 
extern double Haldane(double rr); 
extern double iHaldane(double mm); 
extern double Morgan(double rr); 
extern double iMorgan(double mm); 
extern double CarterFalconer(double rr); 
extern double iCarterFalconer(double mm); 
extern double Rao(double rr); 
extern double iRao(double mm); 
extern double Sturt(double rr); 
extern double iSturt(double mm); 
extern double Felsenstein(double rr); 
extern double iFelsenstein(double mm); 
extern double Karlin(double rr); 
extern double iKarlin(double mm); 
extern int isfile(char *filename);
extern int get_int(void); 
extern int get_next_line(char *buffer, int nn, FILE *fileptr);
extern int get_next_token(char *xtemp, int nn, FILE *fileptr);
extern int dtranspose(double **mm1, double **mm2, int lr, int lc, int ur, int uc);      /*  Transpose a matrix of doubles 
                           **mm1, **mm2, lr,lc,ur,uc     */
extern long get_a_seed(void);
extern char *asctime2(void); 
extern void print_histogram(double *y, int n, int intervals, char *outfile);
extern long get_identifier(char *thefile);
extern void gnuplot_values(double **xx, double **yy, int lr, int ur, int lc, int uc, char *thefile, char *title, char *xaxis, char *yaxis);

#if defined(DSIGN)
extern double dsign(double val1, double val2);        /* v1 = dsign(v1,v2) transfers sign from v2 to v1 */
#endif

#if defined(DIVT)
typedef struct {
  int     quot;
  int     rem;
}       div_t;

typedef struct {
  long     quot;
  long     rem;
}       ldiv_t;
div_t   div(int numer, int denom);
ldiv_t  ldiv(long numer,long denom);
#endif


#if defined(ITOA) 
char   *strupr(char *s);   /* transforms a string to upper case */
char   *strlwr(char *s);   /* transforms a string to lower case */
#endif


extern void pause(void);             /* read characters until a <CR> */
extern void get_field(int xfield, char *xtemp, char *xbuffer);         /* gets the specified field in a string */


extern double ranf(long int inix);
extern double *ran_arry(int size);

extern int      *ivector(int nl, int nh);
extern long     *lvector(int nl, int nh);
extern float    *vector(int nl, int nh);
extern double   *dvector(int nl, int nh);
extern char     *cvector(int nl, int nh);
extern int     **imatrix(int nrl, int nrh, int ncl, int nch);
extern short   **smatrix(int nrl, int nrh, int ncl, int nch);
extern char    **cmatrix(int nrl, int nrh, int ncl, int nch);
extern float   **matrix(int nrl, int nrh, int ncl, int nch);
extern double  **dmatrix(int nrl, int nrh, int ncl, int nch);
extern void    free_ivector(int *v, int nl, int nh);
extern void    free_lvector(long int *v, int nl, int nh);
extern void    free_vector(float *v, int nl, int nh);
extern void    free_dvector(double *v, int nl, int nh);
extern void    free_cvector(char *v, int nl, int nh);
extern void    free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
extern void    free_smatrix(short int **m, int nrl, int nrh, int ncl, int nch);
extern void    free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch);
extern void    free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
extern void    free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
extern void    nrerror(char *error_text);

extern long iran(long int xix, long int in);
extern double gamnl1(double beta, double aa, double bb, double pp, long int xix); 
extern double gamgbh(double beta, long int xix);

#endif

