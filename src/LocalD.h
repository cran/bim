#ifndef LOCALD_GUARD
#define LOCALD_GUARD

/* Local definitions for the compiler */
    /* Define this for Macintosh (CodeWarrior)*/
/*#ifndef THINK_C
#define THINK_C 1
#endif
*/    /*Define this if you need the sign() function*/
#ifndef DSIGN
#define DSIGN 1
#endif
    /*Define this if you need div_t, ldiv_t*/
/*
#ifndef DIVT
#define DIVT 1 
#endif
*/
    /*Define this if you need itoa(), strlwr() and strupr()*/

#ifndef ITOA
#define ITOA 1 
#endif

    /*Define this if using Borland compiler*/
 
/*#ifndef BORLAND
#define BORLAND 1 
#endif*/
 
    /*Define this if using  WinCodewarrior compiler*/
/*#ifndef WINWARRIOR
#define WINWARRIOR 1 
#endif*/

    /* 0 => no debugging output, 1 => memory allocation output */
#ifndef DEBUGGING
#define DEBUGGING 0
#endif

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif
/*Macintosh, Windows and UNIX delimiters*/
/*
#define FILESEP ':'
#define FILESEP '\\'
*/
#define FILESEP '/'

#endif

