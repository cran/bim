#ifndef HELPFILE_GUARD
#define HELPFILE_GUARD
/*
 * $Id: help.h%v 3.50 1993/07/09 05:35:24 woo Exp $
 *
 */


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               Following Modification by C. Basten 6 Sept. 1994
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define HELPFILE "qtlcart.hlp"

#ifndef MAXNAME
#define MAXNAME 64
#endif

#if defined(__EMX__) || defined(DJGPP) || defined(DOS386)
#ifdef MSDOS
#undef MSDOS	/* we have plenty of memory under __EMX__ or DJGPP */
#endif
#ifdef unix
#undef unix	/* we are not unix */
#endif
#endif

#ifdef OS2
  /* GCC defines unix, but no PAGER, so... */
#ifdef unix
#undef unix
#endif
#endif  /* OS2 */

#define	KEYFLAG	'?'	/* leading char in help file topic lines */

#ifndef BORLAND
typedef unsigned char boolean;
#endif


#ifdef	WDLEN
#  define	PATHSIZE	WDLEN
#else
#  define	PATHSIZE	BUFSIZ
#endif


#ifndef TRUE
#define TRUE (1)
#define FALSE (0)
#endif

typedef struct line_s LINEBUF;
struct line_s {
    char *line;			/* the text of this line */
    LINEBUF *next;			/* the next line */
};

typedef struct linkey_s LINKEY;
struct linkey_s {
    char *key;				/* the name of this key */
    long pos;			    /* ftell position */
    LINEBUF *text;			/* the text for this key */
    boolean primary;		/* TRUE -> is a primary name for a text block */
    LINKEY *next;			/* the next key in linked list */
};

typedef struct key_s KEY;
struct key_s {
    char *key;				/* the name of this key */
    long pos;			    /* ftell position */
    LINEBUF *text;			/* the text for this key */
    boolean primary;		/* TRUE -> is a primary name for a text block */
};


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               End of Modification by C. Basten 6 Sept. 1994
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



/* GNUPLOT - help.h */
/*
 * Copyright (C) 1986 - 1993   Thomas Williams, Colin Kelley
 *
 * Permission to use, copy, and distribute this software and its
 * documentation for any purpose with or without fee is hereby granted, 
 * provided that the above copyright notice appear in all copies and 
 * that both that copyright notice and this permission notice appear 
 * in supporting documentation.
 *
 * Permission to modify the software is granted, but not the right to
 * distribute the modified code.  Modifications are to be distributed 
 * as patches to released version.
 *  
 * This software is provided "as is" without express or implied warranty.
 * 
 *
 * AUTHORS
 * 
 *   Original Software:
 *     Thomas Williams,  Colin Kelley.
 * 
 *   Gnuplot 2.0 additions:
 *       Russell Lang, Dave Kotz, John Campbell.
 *
 *   Gnuplot 3.0 additions:
 *       Gershon Elber and many others.
 * 
 * There is a mailing list for gnuplot users. Note, however, that the
 * newsgroup 
 *	comp.graphics.gnuplot 
 * is identical to the mailing list (they
 * both carry the same set of messages). We prefer that you read the
 * messages through that newsgroup, to subscribing to the mailing list.
 * (If you can read that newsgroup, and are already on the mailing list,
 * please send a message info-gnuplot-request@dartmouth.edu, asking to be
 * removed from the mailing list.)
 *
 * The address for mailing to list members is
 *	   info-gnuplot@dartmouth.edu
 * and for mailing administrative requests is 
 *	   info-gnuplot-request@dartmouth.edu
 * The mailing list for bug reports is 
 *	   bug-gnuplot@dartmouth.edu
 * The list of those interested in beta-test versions is
 *	   info-gnuplot-beta@dartmouth.edu
 */

/* Exit status returned by help() */
#define	H_FOUND		0	/* found the keyword */
#define	H_NOTFOUND	1	/* didn't find the keyword */
#define	H_ERROR		(-1)	/* didn't find the help file */


extern void FreeHelp();		/* use this if you need memory */

extern void int_error(char* error_text, int dummy);
extern void get_help(char* helpfile);
extern int help(char* keyword, char* path, boolean* subtopics);
extern int LoadHelp(char*);
extern KEY* FindHelp(char*);
extern void PrintHelp(KEY*, boolean*);
extern LINKEY *storekey(char*);
extern LINEBUF * storeline(char* text);
extern void sortkeys();
extern int instring(char *str, char c);
extern int keycomp(const void* a, const void* b);
extern boolean Ambiguous(KEY* key, int len);
extern void ShowSubtopics(KEY* key, boolean* subtopics);
extern void StartOutput();
extern void OutLine(char*);
extern void EndOutput();
extern KEY *FindHelp(char*);

#endif

