#ifndef PARAMS_GUARD
#define PARAMS_GUARD
/*  
          params.h
          
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


extern void update_opts(char **opt,char  **opt_v,char  **opt_e, int nopts, params *theparams, int flag);
extern void update_params(char **opt_v, int nopts, params *theparams);
extern void set_workdir(params *theparams);

extern void put_hline(FILE *fptr, char ch, int length);
extern void write_trailer(params *theparams,char *chptr,int pfrw);
extern void shift_fn(char *inbuff);
extern void insert_wd(char *buff,char *wd,char *fn);
extern void unset_workdir(params *theparams);
extern void renew_workdir(params *theparams,  char *xtemp,  char **opt, char  **opt_v,  char **opt_e,  int nopts);
extern void renew_resource(params *theparams,  char *xtemp, char   **opt, char  **opt_v, char  **opt_e,  int nopts);
extern int check_params(params *theparams, markermap *themap, aqtl *theqtls, individual *individs, int automatic, int prog);


extern void renew_stem(char *xtemp,params *theparams);
extern void rewrite_param_file(params *theparams,char *filename);
extern void get_param_file(params *theparams,char *qtlrc);
extern params *create_params(params *aparams,int cd,char *qtlrc);

extern int show_opts(FILE *fptr,char *tptr,char  *prog,char  *purpose,char  **opt,char  **opt_v,char  **opt_e,int nopts,int  oflag, params *theparams);
extern void create_opts(char **opt, char **opt_v, char **opt_e, int nopts);
extern void destroy_opts(char **opt, char **opt_v, char **opt_e, int nopts);
extern int parse_cross(params *theparams,  char *xtemp);
extern int     get_cross(char *inputf, params *theparams);
extern int process_arguments(int argc, char **argv, char *thetime, char *purpose,  int nopts,  params *theparams);

extern int get_file_type(char *filename);
extern int file_to_int(char *xtemp);
extern void write_file_type(int ft,  FILE *fptr);
extern void print_head(char *prog,char *filename,char *chptr,int outmode,int filetype,  params *theparams);
extern void check_directory(char *chptr);
extern void quit_banner(char *stri);

extern int whichprogram;


#endif

