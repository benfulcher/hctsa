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
/*Author: Rainer Hegger. Last modified: Sep 7, 2004 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"
#include <math.h>

#define WID_STR "Estimates the average forecast error for a local\n\t\
constant fit as a function of the neighborhood size."


/*number of boxes for the neighbor search algorithm*/
#define NMAX 256

unsigned int nmax=(NMAX-1);
long **box,*list;
unsigned long *found;
double *error;
double **series;

char eps0set=0,eps1set=0,causalset=0,dimset=0;
char *outfile=NULL,stdo=1;
char *column=NULL;
unsigned int dim=1,embed=2,delay=1;
unsigned int verbosity=0xff;
int STEP=1;
double EPS0=1.e-3,EPS1=1.0,EPSF=1.2;
unsigned long LENGTH=ULONG_MAX,exclude=0,CLENGTH=ULONG_MAX,causal;
char *infile=NULL;

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
  fprintf(stderr,"\t-c columns to read [default: 1,...,# of components]\n");
  fprintf(stderr,"\t-m # of components,embedding dimension [default: 1,2]\n");
  fprintf(stderr,"\t-d delay [default: 1]\n");
  fprintf(stderr,"\t-i iterations [default: length]\n");
  fprintf(stderr,"\t-r neighborhood size to start with [default:"
	  " (interval of data)/1000)]\n");
  fprintf(stderr,"\t-R neighborhood size to end with [default:"
	  " interval of data]\n");
  fprintf(stderr,"\t-f factor to increase size [default: 1.2]\n");
  fprintf(stderr,"\t-s steps to forecast [default: 1]\n");
  fprintf(stderr,"\t-C width of causality window [default: steps]\n");
  fprintf(stderr,"\t-o output file name [default: 'datafile.lm']\n");
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
  if ((out=check_option(in,n,'c','s')) != NULL) {
    column=out;
    dimset=1;
  }
  if ((out=check_option(in,n,'m','2')) != NULL)
    sscanf(out,"%u,%u",&dim,&embed);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(in,n,'i','u')) != NULL)
    sscanf(out,"%lu",&CLENGTH);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    eps0set=1;
    sscanf(out,"%lf",&EPS0);
  }
  if ((out=check_option(in,n,'R','f')) != NULL) {
    eps1set=1;
    sscanf(out,"%lf",&EPS1);
  }
  if ((out=check_option(in,n,'f','f')) != NULL)
    sscanf(out,"%lf",&EPSF);
  if ((out=check_option(in,n,'s','u')) != NULL)
    sscanf(out,"%u",&STEP);
  if ((out=check_option(in,n,'C','u')) != NULL) {
    sscanf(out,"%lu",&causal);
    causalset=1;
  }
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void make_fit(long act,unsigned long number)
{
  double *si,cast;
  long i,j;
  
  for (i=0;i<dim;i++) {
    si=series[i];
    cast=si[found[0]+STEP];
    for (j=1;j<number;j++)
      cast += si[found[j]+STEP];
    cast /= (double)number;
    error[i] += sqr(cast-series[i][act+STEP]);
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  unsigned long actfound;
  unsigned long *hfound;
  long pfound,i,j;
  unsigned long clength;
  double interval,min,maxinterval;
  double epsilon;
  double **hser;
  double avfound,*hrms,*hav,sumerror=0.0;
  FILE *file=NULL;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  if (!causalset)
    causal=STEP;

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+4,(size_t)1));
      sprintf(outfile,"%s.lm",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)9,(size_t)1));
      sprintf(outfile,"stdin.lm");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (column == NULL)
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dim,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dim,column,
				      dimset,verbosity);
  maxinterval=0.0;
  for (i=0;i<dim;i++) {
    rescale_data(series[i],LENGTH,&min,&interval);
    if (interval > maxinterval)
      maxinterval=interval;
  }
  interval=maxinterval;

  check_alloc(list=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(found=(unsigned long*)malloc(sizeof(long)*LENGTH));
  check_alloc(hfound=(unsigned long*)malloc(sizeof(long)*LENGTH));
  check_alloc(box=(long**)malloc(sizeof(long*)*NMAX));
  for (i=0;i<NMAX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*NMAX));
  check_alloc(error=(double*)malloc(sizeof(double)*dim));
  check_alloc(hrms=(double*)malloc(sizeof(double)*dim));
  check_alloc(hav=(double*)malloc(sizeof(double)*dim));
  check_alloc(hser=(double**)malloc(sizeof(double*)*dim));
  
  if (eps0set)
    EPS0 /= interval;
  if (eps1set)
    EPS1 /= interval;

  clength=(CLENGTH <= LENGTH) ? CLENGTH-STEP : LENGTH-STEP;

  if (!stdo) {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(file,"#1. size 2. relative forecast error 3. fraction of points\n"
	    "#4. av neighbors found 5. absolute variance of the points\n");
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  for (epsilon=EPS0;epsilon<EPS1*EPSF;epsilon*=EPSF) {
    pfound=0;
    for (i=0;i<dim;i++)
      error[i]=hrms[i]=hav[i]=0.0;
    avfound=0.0;
    make_multi_box(series,box,list,LENGTH-STEP,NMAX,dim,
		   embed,delay,epsilon);
    for (i=(embed-1)*delay;i<clength;i++) {
      for (j=0;j<dim;j++)
	hser[j]=series[j]+i;
      actfound=find_multi_neighbors(series,box,list,hser,LENGTH,
				    NMAX,dim,embed,delay,epsilon,hfound);
      actfound=exclude_interval(actfound,i-causal+1,i+causal+(embed-1)*delay-1,
				hfound,found);
      if (actfound > 2*(dim*embed+1)) {
	make_fit(i,actfound);
	pfound++;
	avfound += (double)(actfound-1);
	for (j=0;j<dim;j++) {
	  hrms[j] += series[j][i+STEP]*series[j][i+STEP];
	  hav[j] += series[j][i+STEP];
	}
      }
    }
    if (pfound > 1) {
      sumerror=0.0;
      for (j=0;j<dim;j++) {
	hav[j] /= pfound;
	hrms[j]=sqrt(fabs(hrms[j]/(pfound-1)-hav[j]*hav[j]*pfound/(pfound-1)));
	error[j]=sqrt(error[j]/pfound)/hrms[j];
	sumerror += error[j];
      }
    }
    if (stdo) {
      if (pfound > 1) {
	fprintf(stdout,"%e %e ",epsilon*interval,sumerror/(double)dim);
	for (j=0;j<dim;j++)
	  fprintf(stdout,"%e ",error[j]);
	fprintf(stdout,"%e %e\n",(double)pfound/(clength-(embed-1)*delay),
		avfound/pfound);
	fflush(stdout);
      }
    }
    else {
      if (pfound > 1) {
	fprintf(file,"%e %e ",epsilon*interval,sumerror/(double)dim);
	for (j=0;j<dim;j++)
	  fprintf(file,"%e ",error[j]);
	fprintf(file,"%e %e\n",(double)pfound/(clength-(embed-1)*delay),
		avfound/pfound);
	fflush(file);
      }
    }
  }
  if (!stdo)
    fclose(file);

  free(list);
  free(hfound);
  free(error);
  free(hrms);
  free(hav);
  free(hser);
  for (i=0;i<NMAX;i++)
    free(box[i]);
  free(box);

  return 0;
}
