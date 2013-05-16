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
/*Author: Rainer Hegger. Last modified: Sep 4, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the average cross forecast error for a zeroth\n\t\
order fit between two series given as two columns of one file."

#ifndef _MATH_H
#include <math.h>
#endif

/*number of boxes for the neighbor search algorithm*/
#define NMAX 128

unsigned int nmax=(NMAX-1);
long **box,*list;
unsigned long *found;
double *series1,*series2;
double interval,min,epsilon;

char epsset=0;
char *infile=NULL;
char *outfile=NULL,stdo=1;
char *COLUMNS=NULL;
unsigned int DIM=3,DELAY=1;
unsigned int verbosity=0xff;
int MINN=30,STEP=1;
double EPS0=1.e-3,EPSF=1.2;
unsigned long LENGTH=ULONG_MAX,exclude=0,CLENGTH=ULONG_MAX;

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
  fprintf(stderr,"\t-x # of lines to be ignored [default: 0]\n");
  fprintf(stderr,"\t-c columns to read [default: 1,2]\n");
  fprintf(stderr,"\t-m embedding dimension [default: 3]\n");
  fprintf(stderr,"\t-d delay [default: 1]\n");
  fprintf(stderr,"\t-n # of reference points [default: length]\n");
  fprintf(stderr,"\t-k minimal number of neighbors for the fit "
	  "[default: 30]\n");
  fprintf(stderr,"\t-r neighborhoud size to start with "
	  "[default: (data interval)/1000]\n");
  fprintf(stderr,"\t-f factor to increase size [default: 1.2]\n");
  fprintf(stderr,"\t-s steps to forecast [default: 1]\n");
  fprintf(stderr,"\t-o output file [default: 'datafile.cze',"
	  " without -o: stdout]\n");
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
    sscanf(out,"%lu",&LENGTH);
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'c','s')) != NULL)
    COLUMNS=out;
  if ((out=check_option(in,n,'m','u')) != NULL)
    sscanf(out,"%u",&DIM);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(in,n,'n','u')) != NULL)
    sscanf(out,"%lu",&CLENGTH);
  if ((out=check_option(in,n,'k','u')) != NULL)
    sscanf(out,"%u",&MINN);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&EPS0);
  }
  if ((out=check_option(in,n,'f','f')) != NULL)
    sscanf(out,"%lf",&EPSF);
  if ((out=check_option(in,n,'s','u')) != NULL)
    sscanf(out,"%u",&STEP);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

double make_fit(unsigned long act,unsigned long number,unsigned long istep)
{
  double casted=0.0;
  int i;
  
  for (i=0;i<number;i++)
    casted += series1[found[i]+istep];
  casted /= number;

  return (casted-series2[act+istep])*(casted-series2[act+istep]);
}

int main(int argc,char **argv)
{
  char stdi=0;
  char alldone,*done;
  unsigned long i,j,actfound;
  unsigned long clength;
  unsigned int dummy=2;
  double rms2,av2,*error;
  double **both,hinter;
  FILE *file;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,0L,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      sprintf(outfile,"%s.cze",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.cze");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (COLUMNS == NULL)
    both=(double**)get_multi_series(infile,&LENGTH,exclude,&dummy,"",
				    (char)1,verbosity);
  else
    both=(double**)get_multi_series(infile,&LENGTH,exclude,&dummy,COLUMNS,
				    (char)0,verbosity);
  series1=both[0];
  series2=both[1];
  rescale_data(series1,LENGTH,&min,&hinter);
  interval=hinter;
  rescale_data(series2,LENGTH,&min,&hinter);
  interval=(interval+hinter)/2.0;

  variance(series2,LENGTH,&av2,&rms2);
  
  check_alloc(list=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(found=(unsigned long*)malloc(sizeof(long)*LENGTH));
  check_alloc(done=(char*)malloc(sizeof(char)*LENGTH));
  check_alloc(box=(long**)malloc(sizeof(long*)*NMAX));
  check_alloc(error=(double*)malloc(sizeof(double)*STEP));
  for (i=0;i<STEP;i++)
    error[i]=0.0;

  for (i=0;i<NMAX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*NMAX));
    
  for (i=0;i<LENGTH;i++)
    done[i]=0;

  alldone=0;
  if (epsset)
    EPS0 /= interval;

  epsilon=EPS0/EPSF;
  clength=(CLENGTH <= LENGTH) ? CLENGTH-STEP : LENGTH-STEP;

  while (!alldone) {
    alldone=1;
    epsilon*=EPSF;
    make_box(series1,box,list,LENGTH-STEP,NMAX,DIM,DELAY,epsilon);
    for (i=(DIM-1)*DELAY;i<clength;i++)
      if (!done[i]) {
	actfound=find_neighbors(series1,box,list,series2+i,LENGTH,NMAX,
				DIM,DELAY,epsilon,found);
	if (actfound >= MINN) {
	  for (j=1;j<=STEP;j++)
	    error[j-1] += make_fit(i,actfound,j);
	  done[i]=1;
	}
	alldone &= done[i];
      }
  }
  if (stdo) {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    for (i=0;i<STEP;i++)
      fprintf(stdout,"%lu %e\n",i+1,
	      sqrt(error[i]/(clength-(DIM-1)*DELAY))/rms2);
  }
  else {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    for (i=0;i<STEP;i++)
      fprintf(file,"%lu %e\n",i+1,
	      sqrt(error[i]/(clength-(DIM-1)*DELAY))/rms2);
    fclose(file);
  }
  
  return 0;
}
