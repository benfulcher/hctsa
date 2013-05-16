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
/*Author: Rainer Hegger. Last modified: Aug 27, 2004 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the average forecast error for a zeroth\n\t\
order fit from a multidimensional time series"


#ifndef _MATH_H
#include <math.h>
#endif

/*number of boxes for the neighbor search algorithm*/
#define NMAX 512

unsigned int nmax=(NMAX-1);
long **box,*list;
unsigned long *found;
double **series,**diffs;
double interval,min,epsilon;

char epsset=0,dimset=0,clengthset=0,causalset=0;
char *infile=NULL;
char *outfile=NULL,stdo=1;
char *COLUMNS=NULL;
unsigned int embed=2,dim=1,DELAY=1,MINN=30;
unsigned long STEP=1,causal;
unsigned int verbosity=0x1;
double EPS0=1.e-3,EPSF=1.2;
unsigned long refstep=1;
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
  fprintf(stderr,"\t-c columns to read [default: 1,...,X]\n");
  fprintf(stderr,"\t-m dimension and embedding dimension"
	  " [default: %d,%d]\n",dim,embed);
  fprintf(stderr,"\t-d delay [default: %d]\n",DELAY);
  fprintf(stderr,"\t-n # of reference points [default: length]\n");
  fprintf(stderr,"\t-S temporal distance between the reference points"
	  " [default: %lu]\n",refstep);
  fprintf(stderr,"\t-k minimal number of neighbors for the fit "
	  "[default: %d]\n",MINN);
  fprintf(stderr,"\t-r neighborhoud size to start with "
	  "[default: (data interval)/1000]\n");
  fprintf(stderr,"\t-f factor to increase size [default: 1.2]\n");
  fprintf(stderr,"\t-s steps to forecast [default: 1]\n");
  fprintf(stderr,"\t-C width of causality window [default: steps]\n");
  fprintf(stderr,"\t-o output file [default: 'datafile.zer',"
	  " without -o: stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
	  "2='give individual forecast errors for the max step'\n");
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
  if ((out=check_option(in,n,'m','2')) != NULL) {
    dimset=1;
    sscanf(out,"%u%*c%u",&dim,&embed);
    if (embed == 0)
      embed=1;
  }
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(in,n,'n','u')) != NULL) {
    sscanf(out,"%lu",&CLENGTH);
    clengthset=1;
  }
  if ((out=check_option(in,n,'S','u')) != NULL)
    sscanf(out,"%lu",&refstep);
  if ((out=check_option(in,n,'k','u')) != NULL)
    sscanf(out,"%u",&MINN);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&EPS0);
  }
  if ((out=check_option(in,n,'f','f')) != NULL)
    sscanf(out,"%lf",&EPSF);
  if ((out=check_option(in,n,'s','u')) != NULL)
    sscanf(out,"%lu",&STEP);
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

void make_fit(long act,unsigned long number,long istep,double **error)
{
  double casted,*help;
  long i,j,h;
  
  h=istep-1;
  for (j=0;j<dim;j++) {
    casted=0.0;
    help=series[j]+istep;
    for (i=0;i<number;i++)
      casted += help[found[i]];
    casted /= number;
    diffs[j][act]=casted-help[act];
    error[j][h] += sqr(casted-help[act]);
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  char alldone,*done;
  long i,j,hi;
  unsigned long *hfound;
  unsigned long actfound;
  unsigned long clength;
  double *rms,*av,**error,**hser,*hinter;
  FILE *file;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);

  if ((2*STEP+causal) >= ((long)LENGTH-(long)(embed*DELAY)-(long)MINN)) {
    fprintf(stderr,"steps to forecast (-s) too large. Exiting!\n");
    exit(ZEROTH__STEP_TOO_LARGE);
  }
  if (!causalset)
    causal=STEP;

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
      sprintf(outfile,"%s.zer",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.zer");
    }
  }
  if (!stdo)
    test_outfile(outfile);
  
  if (COLUMNS == NULL)
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dim,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dim,COLUMNS,
				      dimset,verbosity);

  check_alloc(hser=(double**)malloc(sizeof(double*)*dim));
  check_alloc(av=(double*)malloc(sizeof(double)*dim));
  check_alloc(rms=(double*)malloc(sizeof(double)*dim));
  check_alloc(hinter=(double*)malloc(sizeof(double)*dim));
  interval=0.0;
  for (i=0;i<dim;i++) {
    rescale_data(series[i],LENGTH,&min,&hinter[i]);
    variance(series[i],LENGTH,&av[i],&rms[i]);
    interval += hinter[i];
  }
  interval /= (double)dim;

  check_alloc(list=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(found=(unsigned long*)malloc(sizeof(long)*LENGTH));
  check_alloc(hfound=(unsigned long*)malloc(sizeof(long)*LENGTH));
  check_alloc(done=(char*)malloc(sizeof(char)*LENGTH));
  check_alloc(box=(long**)malloc(sizeof(long*)*NMAX));
  check_alloc(error=(double**)malloc(sizeof(double*)*dim));
  check_alloc(diffs=(double**)malloc(sizeof(double*)*dim));
  for (j=0;j<dim;j++) {
    check_alloc(diffs[j]=(double*)malloc(sizeof(double)*LENGTH));
    check_alloc(error[j]=(double*)malloc(sizeof(double)*STEP));
    for (i=0;i<STEP;i++)
      error[j][i]=0.0;
  }
  
  for (i=0;i<NMAX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*NMAX));
    
  for (i=0;i<LENGTH;i++)
    done[i]=0;

  alldone=0;
  if (epsset)
    EPS0 /= interval;

  epsilon=EPS0/EPSF;

  if (!clengthset)
    CLENGTH=LENGTH;
  clength=((CLENGTH*refstep+STEP) <= LENGTH) ? CLENGTH : 
    (LENGTH-(long)STEP)/refstep;

  while (!alldone) {
    alldone=1;
    epsilon*=EPSF;
    make_multi_box(series,box,list,LENGTH-(long)STEP,NMAX,(unsigned int)dim,
		   (unsigned int)embed,(unsigned int)DELAY,epsilon);
    for (i=(embed-1)*DELAY;i<clength;i++)
      if (!done[i]) {
	hi=i*refstep;
	for (j=0;j<dim;j++)
	  hser[j]=series[j]+hi;
	actfound=find_multi_neighbors(series,box,list,hser,LENGTH,NMAX,
				       (unsigned int)dim,(unsigned int)embed,
				       (unsigned int)DELAY,epsilon,hfound);
	actfound=exclude_interval(actfound,hi-(long)causal+1,
				  hi+causal+(embed-1)*DELAY-1,hfound,found);	
	if (actfound >= MINN) {
	  for (j=1;j<=STEP;j++) {
	    make_fit(hi,actfound,j,error);
	  }
	  done[i]=1;
	}
	alldone &= done[i];
      }
  }
  if (stdo) {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    for (i=0;i<STEP;i++) {
      if (verbosity&VER_USR1)
	fprintf(stdout,"# %lu ",i+1);
      else
	fprintf(stdout,"%lu ",i+1);
      for (j=0;j<dim;j++) 
	fprintf(stdout,"%e ",
		sqrt(error[j][i]/(clength-(embed-1)*DELAY))/rms[j]);
      fprintf(stdout,"\n");
    }
    if (verbosity&VER_USR1) {
      for (i=(embed-1)*DELAY;i<clength;i++) {
	hi=i*refstep;
	for (j=0;j<dim;j++)
	  fprintf(stdout,"%e ",diffs[j][hi]*hinter[j]);
	fprintf(stdout,"\n");
      }
    }
  }
  else {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    for (i=0;i<STEP;i++) {
      if (verbosity&VER_USR1)
	fprintf(file,"# %lu ",i+1);
      else
	fprintf(file,"%lu ",i+1);
      for (j=0;j<dim;j++) 
	fprintf(file,"%e ",sqrt(error[j][i]/(clength-(embed-1)*DELAY))/rms[j]);
      fprintf(file,"\n");
    }
    if (verbosity&VER_USR1) {
      for (i=(embed-1)*DELAY;i<clength;i++) {
	hi=i*refstep;
	for (j=0;j<dim;j++)
	  fprintf(file,"%e ",diffs[j][hi]*hinter[j]);
	fprintf(file,"\n");
      }
    }
    fclose(file);
  }

  if (outfile != NULL)
    free(outfile);
  if (infile != NULL)
    free(infile);
  if (COLUMNS != NULL)
    free(COLUMNS);
  for (i=0;i<dim;i++) {
    free(series[i]);
    free(diffs[i]);
    free(error[i]);
  }
  free(series);
  free(diffs);
  free(hser);
  free(error);
  free(av);
  free(rms);
  free(list);
  free(found);
  free(hfound);
  free(done);
  for (i=0;i<NMAX;i++)
    free(box[i]);
  free(box);

  return 0;
}
