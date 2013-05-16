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
/* Author: Rainer Hegger Last modified: Feb 22, 2006 */
/* Changes: 
   02/22/06: Remove this strange else in start_box that 
             did not compile anyways
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the Renyi entropy of Qth order\n\t\
using a partition instead of a covering."

typedef struct {
  double *hist;
  void *ptr;
} hliste;

unsigned long LENGTH=ULONG_MAX,exclude=0;
unsigned int maxembed=10,dimension=1,DELAY=1,EPSCOUNT=20;
unsigned int verbosity=0xff;
double Q=2.0,EPSMIN=1.e-3,EPSMAX=1.0;
char dimset=0,epsminset=0,epsmaxset=0;
char *outfile=NULL;
char *column=NULL;

int epsi;
unsigned long length;
double EPSFAKTOR;
unsigned int **which_dims;
double *histo;
double **series;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"Usage: %s [Options]\n",progname);
  fprintf(stderr,"Options:\n");
  fprintf(stderr,"\t-l # of datapoints [Default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [Default: %lu]\n",exclude);
  fprintf(stderr,"\t-M # of columns,maximal embedding dimension "
	  "[Default: %u,%u]\n",dimension,maxembed);
  fprintf(stderr,"\t-c columns to read  [Default: 1,...,#of compon.]\n");
  fprintf(stderr,"\t-d delay [Default: %u]\n",DELAY);
  fprintf(stderr,"\t-Q order of the Renyi entropy [Default: %.1f]\n",Q);
  fprintf(stderr,"\t-r minimal epsilon [Default: (data interval)/1000]\n");
  fprintf(stderr,"\t-R maximal epsilon [Default: data interval]\n");
  fprintf(stderr,"\t-# # of epsilons to use [Default: %u]\n",EPSCOUNT);
  fprintf(stderr,"\t-o output file name [Default: 'datafile'.box]\n");
  fprintf(stderr,"\t-V verbosity level [Default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n\n");
  exit(0);
}

void scan_options(int n,char **in)
{
  char *out;
  
  if ((out=check_option(in,n,'l','u')) != NULL)
    sscanf(out,"%lu",&LENGTH);
  if ((out=check_option(in,n,'c','s')) != NULL)
    column=out;
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'M','2')) != NULL) {
    sscanf(out,"%u,%u",&dimension,&maxembed);
    dimset=1;
  }
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(in,n,'Q','f')) != NULL)
    sscanf(out,"%lf",&Q);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    sscanf(out,"%lf",&EPSMIN);
    epsminset=1;
  }
  if ((out=check_option(in,n,'R','f')) != NULL) {
    sscanf(out,"%lf",&EPSMAX);
    epsmaxset=1;
  }
  if ((out=check_option(in,n,'#','u')) != NULL)
    sscanf(out,"%u",&EPSCOUNT);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'o','s')) != NULL)
    outfile=out;
}

hliste *make_histo(void)
{
  int i;
  hliste *element;
  
  check_alloc(element=(hliste*)malloc(sizeof(hliste)));
  element->ptr=NULL;
  check_alloc(element->hist=(double*)malloc(sizeof(double)*maxembed*dimension));
  for (i=0;i<maxembed*dimension;i++)
    element->hist[i]=0.0;
  
  return element;
}

void next_dim(int wd,int n,unsigned int *first)
{
  int i,which,d1,comp;
  double epsinv,norm,p;
  unsigned int **act;
  int *found,hf;

  comp=which_dims[wd][0];
  d1=which_dims[wd][1]*DELAY;

  epsinv=(double)epsi;
  norm=(double)length;

  check_alloc(act=(unsigned int**)malloc(epsi*sizeof(int*)));
  check_alloc(found=(int*)malloc(epsi*sizeof(int)));
  
  for (i=0;i<epsi;i++) {
    found[i]=0;
    act[i]=NULL;
  }
  
  for (i=0;i<n;i++) {
    which=(int)(series[comp][first[i]+d1]*epsinv);
    hf= ++found[which];
    check_alloc(act[which]=
		realloc((unsigned int*)act[which],hf*sizeof(unsigned int)));
    act[which][hf-1]=first[i];
  }
  
  for (i=0;i<epsi;i++)
    if (found[i]) {
      p=(double)(found[i])/(norm);
      if (Q == 1.0)
	histo[wd] -= p*log(p);
      else
	histo[wd] += pow(p,Q);
    }
  
  if (wd<(maxembed*dimension-1))
    for (i=0;i<epsi;i++)
      if (found[i])
	next_dim(wd+1,found[i],act[i]);
  
  for (i=0;i<epsi;i++)
    if (found[i])
      free(act[i]);
  
  free(act);
  free(found);
}

void start_box(void)
{
  int i,which;
  double epsinv,norm,p;
  unsigned int **act;
  int *found,hf;
  void next_dim();
  
  epsinv=(double)epsi;
  norm=(double)length;
  
  check_alloc(act=(unsigned int**)malloc(epsi*sizeof(int*)));
  check_alloc(found=(int*)malloc(epsi*sizeof(int)));
  
  for (i=0;i<epsi;i++) {
    found[i]=0;
    act[i]=NULL;
  }
  
  for (i=0;i<length;i++) {
    which=(int)(series[0][i]*epsinv);
    hf= ++found[which];
    check_alloc(act[which]=
		realloc((unsigned int*)act[which],hf*sizeof(unsigned int)));
    act[which][hf-1]=i;
  }
  
  for (i=0;i<epsi;i++)
    if (found[i]) {
      p=(double)(found[i])/(norm);
      if (Q == 1.0)
	histo[0] -= p*log(p);
      else
	histo[0] += pow(p,Q);
    }
  
  if (1<dimension*maxembed) {
    for (i=0;i<epsi;i++) {
      if (found[i])
	next_dim(1,found[i],act[i]);
    }
  }
  /*
  else {
    if (1<maxembed)
      for (i=0;i<epsi;i++) {
	if (found[i])
	  next_dim(1,found[i],act[i]);
      }
  }
  */

  for (i=0;i<epsi;i++)
    if (found[i])
      free(act[i]);
  
  free(act);
  free(found);
}

int main(int argc,char **argv)
{
  int i,j,k,count,epsi_old=0,epsi_test;
  void *root;
  hliste *histo_el;
  double *deps,heps;
  double min,interval,maxinterval;
  char *infile=NULL,stdi=0;
  FILE *fHq;


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
      sprintf(outfile,"%s.box",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.box");
    }
  }
  test_outfile(outfile);

  if (column == NULL)
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dimension,"",
				      dimset,verbosity);
  else
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dimension,
				      column,dimset,verbosity);
  maxinterval=0.0;
  for (i=0;i<dimension;i++) {
    rescale_data(series[i],LENGTH,&min,&interval);
    if (interval > maxinterval)
      maxinterval=interval;
  }
  if (epsminset)
    EPSMIN /= maxinterval;
  if (epsmaxset)
    EPSMAX /= maxinterval;
  for (i=0;i<dimension;i++) {
    for (j=0;j<LENGTH;j++)
      if (series[i][j] >= 1.0)
	series[i][j] -= EPSMIN/2.0;
  }

  check_alloc(histo=(double*)malloc(sizeof(double)*maxembed*dimension));
  check_alloc(deps=(double*)malloc(sizeof(double)*EPSCOUNT));
  check_alloc(which_dims=(unsigned int**)malloc(sizeof(int*)*
						maxembed*dimension));
  for (i=0;i<maxembed*dimension;i++)
    check_alloc(which_dims[i]=(unsigned int*)malloc(sizeof(int)*2));
  for (i=0;i<maxembed;i++)
    for (j=0;j<dimension;j++) {
      which_dims[i*dimension+j][0]=j;
      which_dims[i*dimension+j][1]=i;
    }
  
  histo_el=make_histo();
  root=histo_el;
  
  if (EPSCOUNT >1)
    EPSFAKTOR=pow(EPSMAX/EPSMIN,1.0/(double)(EPSCOUNT-1));
  else
    EPSFAKTOR=1.0;

  length=LENGTH-(maxembed-1)*DELAY;

  heps=EPSMAX*EPSFAKTOR;
  
  for (k=0;k<EPSCOUNT;k++) {
    count++;
    for (i=0;i<maxembed*dimension;i++)
      histo[i]=0.0;
    do {
      heps /= EPSFAKTOR;
      epsi_test=(int)(1./heps);
    } while (epsi_test <= epsi_old);
    
    epsi=epsi_test;
    epsi_old=epsi;
    deps[k]=heps;
    
    start_box();
    histo_el=root;
    while (histo_el->ptr != NULL)
      histo_el=histo_el->ptr;
    
    for (i=0;i<maxembed*dimension;i++)
      if (Q == 1.0)
	histo_el->hist[i]=histo[i];
      else
	histo_el->hist[i]=log(histo[i])/(1.0-Q);
    
    histo_el->ptr=make_histo();
    histo_el=histo_el->ptr;
    fHq=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);

    for (i=0;i<maxembed*dimension;i++) {
      fprintf(fHq,"#component = %d embedding = %d\n",which_dims[i][0]+1,
	      which_dims[i][1]+1);
      histo_el=root;
      for (j=0;j<=k;j++) {
	if (i == 0)
	  fprintf(fHq,"%e %e %e\n",deps[j]*maxinterval,
		  histo_el->hist[i],histo_el->hist[i]);
	else
	  fprintf(fHq,"%e %e %e\n",deps[j]*maxinterval,
		  histo_el->hist[i],histo_el->hist[i]-histo_el->hist[i-1]);
	histo_el=histo_el->ptr;
      }
      fprintf(fHq,"\n");
    }
    fclose(fHq);
  }

  return 0;
}
