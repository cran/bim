/* Prototypes for all user accessible RANLIB routines */
#define FLOAT double

extern void advnst(long k);
extern FLOAT genbet(FLOAT aa,FLOAT bb);
extern FLOAT genchi(FLOAT df);
extern FLOAT genexp(FLOAT av);
extern FLOAT genf(FLOAT dfn, FLOAT dfd);
extern FLOAT gengam(FLOAT a,FLOAT r);
extern void genmn(FLOAT *parm,FLOAT *x,FLOAT *work);
extern void genmul(long n,FLOAT *p,long ncat,long *ix);
extern void genmul2(long n,long ncat,long *ix);
extern FLOAT gennch(FLOAT df,FLOAT xnonc);
extern FLOAT gennf(FLOAT dfn, FLOAT dfd, FLOAT xnonc);
extern FLOAT gennor(FLOAT av,FLOAT sd);
extern void genprm(long *iarray,int larray);
extern FLOAT genunf(FLOAT low,FLOAT high);
extern void getsd(long *iseed1,long *iseed2);
extern void gscgn(long getset,long *g);
extern long ignbin(long n,FLOAT pp);
extern long ignnbn(long n,FLOAT p);
extern long ignlgi(void);
extern long ignpoi(FLOAT mu);
extern long ignuin(long low,long high);
extern void initgn(long isdtyp);
extern long mltmod(long a,long s,long m);
extern void phrtsd(char* phrase,long* seed1,long* seed2);
extern FLOAT RANF(void);
extern void setall(long iseed1,long iseed2);
extern void setant(long qvalue);
extern void setgmn(FLOAT *meanv,FLOAT *covm,long p,FLOAT *parm);
extern void setsd(long iseed1,long iseed2);
extern FLOAT sexpo(void);
extern FLOAT sgamma(FLOAT a);
extern FLOAT snorm(void);
