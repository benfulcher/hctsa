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
/*Author: Rainer Hegger, Last modified: Mar 20, 1999 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Make a Poincare section"

char *outfile=NULL,dimset=0,compset=0,whereset=0,stdo=1;
char *infile=NULL;
unsigned long length=ULONG_MAX,count,exclude=0;
int dim=2,comp=2,delay=1,dir=0;
unsigned int column=1;
unsigned int verbosity=0xff;
double *series,min,max,average=0.0,where;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [Options]\n\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of points to be used [Default: whole file]\n");
  fprintf(stderr,"\t-x #of lines to be ignored [Default: 0]\n");
  fprintf(stderr,"\t-c column to read  [Default: 1]\n");
  fprintf(stderr,"\t-m embedding dimension [Default: 2]\n");
  fprintf(stderr,"\t-d delay [Default: 1]\n");
  fprintf(stderr,"\t-q component to cut [Default: last]\n");
  fprintf(stderr,"\t-C direction of the cut (0: from below,1: from above)"
	  "\n\t\t[Default: 0]\n");
  fprintf(stderr,"\t-a set crossing at [Default: average of data]\n");
  fprintf(stderr,"\t-o outfile [Default: 'datafile'.poin]\n");
  fprintf(stderr,"\t-V verbosity level [Default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h  show these options\n");
  fprintf(stderr,"\n");
  exit(0);
}


void scan_options(int n,char** in)
{
  char *out;
  
  if ((out=check_option(in,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(in,n,'m','u')) != NULL) {
    dimset=1;
    sscanf(out,"%u",&dim);
  }
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(in,n,'q','u')) != NULL) {
    compset=1;
    sscanf(out,"%u",&comp);
  }
  if ((out=check_option(in,n,'C','u')) != NULL)
    sscanf(out,"%u",&dir);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'a','f')) != NULL) {
    whereset=1;
    sscanf(out,"%lf",&where);
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void poincare(void)
{
  unsigned long i;
  long j,jd;
  double delta,xcut;
  double time=0.0,lasttime=0.0;
  FILE *fout=NULL;

  if (!stdo) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  if (dir == 0) {
    for (i=(comp-1)*delay;i<length-(dim-comp)*delay-1;i++) {
      if ((series[i] < where) && (series[i+1] >= where)) {
	delta=(series[i]-where)/(series[i]-series[i+1]);
	time=(double)i+delta;
	if (lasttime > 0.0) {
	  for (j= -(comp-1);j<=dim-comp;j++) {
	    if (j != 0) {
	      jd=i+j*delay;
	      xcut=series[jd]+delta*(series[jd+1]-series[jd]);
	      if (!stdo)
		fprintf(fout,"%e ",xcut);
	      else
		fprintf(stdout,"%e ",xcut);
	    }
	  }
	  if (!stdo) {
	    fprintf(fout,"%e\n",time-lasttime);
	    fflush(fout);
	  }
	  else {
	    fprintf(stdout,"%e\n",time-lasttime);
	    fflush(stdout);
	  }
	  count++;
	}
	lasttime=time;
      }
    }
  }
  else {
    for (i=(comp-1)*delay;i<length-(dim-comp)*delay-1;i++) {
      if ((series[i] > where) && (series[i+1] <= where)) {
	delta=(series[i]-where)/(series[i]-series[i+1]);
	time=(double)i+delta;
	if (lasttime > 0.0) {
	  for (j= -(comp-1);j<=dim-comp;j++) {
	    if (j != 0) {
	      jd=i+j*delay;
	      xcut=series[jd]+delta*(series[jd+1]-series[jd]);
	      if (!stdo)
		fprintf(fout,"%e ",xcut);
	      else
		fprintf(stdout,"%e ",xcut);
	    }
	  }
	  if (!stdo) {
	    fprintf(fout,"%e\n",time-lasttime);
	    fflush(fout);
	  }
	  else {
	    fprintf(stdout,"%e\n",time-lasttime);
	    fflush(stdout);
	  }
	  count++;
	}
	lasttime=time;
      }
    }
  }
  if (!stdo)
    fclose(fout);
}

int main(int argc,char** argv)
{
  char stdi=0;
  long i;
  double var;

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
      check_alloc(outfile=(char*)calloc(strlen(infile)+6,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".poin");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)11,(size_t)1));
      strcpy(outfile,"stdin.poin");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  variance(series,length,&average,&var);
  min=max=series[0];
  for (i=1;i<length;i++) {
    if (series[i] < min) min=series[i];
    if (series[i] > max) max=series[i];
  }

  if (!whereset)
    where=average;
  if (dimset && !compset)
    comp=dim;
  
  if (comp > dim) {
    fprintf(stderr,"Component to cut is larger than dimension. Exiting!\n");
    exit(POINCARE__WRONG_COMPONENT);
  }
  if ((where < min) || (where > max)) {
    fprintf(stderr,"You want to cut outside the data interval which is [%e,"
	    "%e]\n",min,max);
    exit(POINCARE__OUTSIDE_REGION);
  }
  poincare();

  return 0;
}
