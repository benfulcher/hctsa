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
/*Author: Rainer Hegger. Last modified: Nov 23, 2000 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "routines/tsa.h"

#define WID_STR "Rescales the data"

unsigned long length=ULONG_MAX,exclude=0;
unsigned int dim=1;
unsigned int verbosity=0xff;
char *column=NULL;
char *outfile=NULL,stdo=1,set_av=0,set_var=0,dimset=0;
char *infile=NULL;
double **series;
double xmin=0.0,xmax=1.0;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to use [default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default: 0]\n");
  fprintf(stderr,"\t-m # of components to be read [default: %u]\n",dim);
  fprintf(stderr,"\t-c columns to read [default: 1,...,# of components]\n");
  fprintf(stderr,"\t-z minimum of the new series [default: 0.0]\n");
  fprintf(stderr,"\t-Z maximum of the new series [default: 1.0]\n");
  fprintf(stderr,"\t-a create a series with average value equals 0\n");
  fprintf(stderr,"\t-v create a series with variance 1\n");
  fprintf(stderr,"\t-o output file name [default: 'datafile'.res']\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **in)
{
  char *out;

  if ((out=check_option(in,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'m','u')) != NULL) {
    sscanf(out,"%u",&dim);
    dimset=1;
  }
  if ((out=check_option(in,n,'c','s')) != NULL)
    column=out;
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'z','f')) != NULL)
    sscanf(out,"%lf",&xmin);
  if ((out=check_option(in,n,'Z','f')) != NULL)
    sscanf(out,"%lf",&xmax);
  if ((out=check_option(in,n,'a','n')) != NULL)
    set_av=1;
  if ((out=check_option(in,n,'v','n')) != NULL)
    set_var=1;
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  FILE *file;
  double min,max;
  double av,varianz;
  long i,n;
    
  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".res");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.res");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (xmin >= xmax) {
    fprintf(stderr,"Choosing the minimum larger or equal the maximum\n"
	    "makes no sense. Exiting!\n");
    exit(RESCALE__WRONG_INTERVAL);
  }

  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,column,
				      dimset,verbosity);

  for (n=0;n<dim;n++) {
    variance(series[n],length,&av,&varianz);
    
    if (set_av)
      for (i=0;i<length;i++)
	series[n][i] -= av;

    if (set_var)
      for (i=0;i<length;i++)
	series[n][i] /= varianz;
  
    if (!set_var && !set_av) {
      rescale_data(series[n],length,&min,&max);
      for (i=0;i<length;i++)
	series[n][i]=series[n][i]*(xmax-xmin)+xmin;
    }
  }

  if (stdo) {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    for (i=0;i<length;i++) {
      fprintf(stdout,"%e",series[0][i]);
      for (n=1;n<dim;n++)
	fprintf(stdout," %e",series[n][i]);
      fprintf(stdout,"\n");
    }
  }
  else {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    for (i=0;i<length;i++) {
      fprintf(file,"%e",series[0][i]);
      for (n=1;n<dim;n++)
	fprintf(file," %e",series[n][i]);
      fprintf(file,"\n");
    }
    fclose(file);
  }

  return 0;
}
