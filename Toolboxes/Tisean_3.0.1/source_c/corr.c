/*
 *   This file is part of TISEAN
 *
 *   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
 *
 *   TISEAN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   TISEAN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TISEAN; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/*Author: Rainer Hegger. Last modified: Sep 3, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the autocorrelations of a data set"

char *format,*outfile=NULL,stout=1,normalize=1;
unsigned int column=1;
unsigned int verbosity=0xff;
unsigned long tau=100,length=ULONG_MAX,exclude=0;
double *array;
double av,var;
char *infile=NULL;

void show_options(char *progname) 
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [Options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l length [default is whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default 0]\n");
  fprintf(stderr,"\t-c column to read [default is 1]\n");
  fprintf(stderr,"\t-D corrlength  [default is 100]\n");
  fprintf(stderr,"\t-n don\'t normalize to the variance"
	  " of the data [not set]\n");
  fprintf(stderr,"\t-o output_file  [default is 'datafile'.cor; no -o"
  " means stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default is 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  fprintf(stderr,"\n");
  exit(0);
}

void scan_options(int argc,char **argv)
{
  char *out;

  if ((out=check_option(argv,argc,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(argv,argc,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(argv,argc,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(argv,argc,'D','u')) != NULL)
    sscanf(out,"%lu",&tau);
  if ((out=check_option(argv,argc,'n','n')) != NULL)
    normalize=0;
  if ((out=check_option(argv,argc,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,argc,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

double corr(long i)
{
  long j;
  double c=0.0;
  
  for (j=0;j<(length-i);j++)
    c += array[j]*array[j+i];

  return c/(length-i);
}

int main(int argc,char** argv)
{
  char stdi=0;
  long i;
  FILE *fout=NULL;
  
  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,&column,verbosity);
  if (infile == NULL)
    stdi=1;
  
  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".cor");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.cor");
    }
  }
  if (!stout)
    test_outfile(outfile);

  array=(double*)get_series(infile,&length,exclude,column,verbosity);

  if (tau >= length)
    tau=length-1;

  variance(array,length,&av,&var);

  if (normalize) {
    for (i=0;i<length;i++)
      array[i] -= av;
  }

  if (!stout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(fout,"# average=%e\n",av);
    fprintf(fout,"# standard deviation=%e\n",var);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"# average=%e\n",av);
    fprintf(stdout,"# standard deviation=%e\n",var);
  }
  if (normalize)
    var *= var;
  else
    var=1.0;

  for (i=0;i<=tau;i++)
    if (!stout) {
      fprintf(fout,"%ld %e\n",i,corr(i)/var);
      fflush(fout);
    }
    else {
      fprintf(stdout,"%ld %e\n",i,corr(i)/var);
      fflush(stdout);
    }
  if (!stout)
    fclose(fout);

  if (outfile != NULL)
    free(outfile);
  if (infile != NULL)
    free(infile);
  free(array);

  return 0;
}
