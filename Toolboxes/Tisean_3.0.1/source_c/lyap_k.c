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
/*Author: Rainer Hegger. Last modified: Sep 3, 1999*/
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the maximal Lyapunov exponent using the Kantz\n\t\
algorithm"

#define BOX 128
const unsigned int ibox=BOX-1;

unsigned long length=ULONG_MAX;
unsigned long exclude=0;
unsigned long reference=ULONG_MAX;
unsigned int maxdim=2;
unsigned int mindim=2;
unsigned int delay=1;
unsigned int column=1;
unsigned int epscount=5;
unsigned int maxiter=50;
unsigned int window=0;
unsigned int verbosity=0xff;
double epsmin=1.e-3,epsmax=1.e-2;
char eps0set=0,eps1set=0;
char *outfile=NULL;
char *infile=NULL;

double *series,**lyap;
long box[BOX][BOX],*liste,**lfound,*found,**count;
double max,min;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);

  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be "
	  "interpreted as a possible datafile.\nIf no datafile "
	  "is given stdin is read. Just - also means stdin\n");
  fprintf(stderr,"\t-l # of data [default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default: 0]\n");
  fprintf(stderr,"\t-c column to read [default: 1]\n");
  fprintf(stderr,"\t-M maxdim [default: 2]\n");
  fprintf(stderr,"\t-m mindim [default: 2]\n");
  fprintf(stderr,"\t-d delay [default: 1]\n");
  fprintf(stderr,"\t-r mineps [default: (data interval)/1000]\n");
  fprintf(stderr,"\t-R maxeps [default: (data interval)/100]\n");
  fprintf(stderr,"\t-# # of eps [default: 5]\n");
  fprintf(stderr,"\t-n # of reference points [default: # of data]\n");
  fprintf(stderr,"\t-s # of iterations [default: 50]\n");
  fprintf(stderr,"\t-t time window [default: 0]\n");
  fprintf(stderr,"\t-o outfile [default: 'datafile'.lyap]\n");
  fprintf(stderr,"\t-V verbosity level [default: 3]\n\t\t"
	  "0='only panic messages'\n\t\t"
	  "1='+ input/output messages'\n\t\t"
	  "2='+ plus statistics'\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **str)
{
  char *out;
  
  if ((out=check_option(str,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(str,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(str,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(str,n,'M','u')) != NULL)
    sscanf(out,"%u",&maxdim);
  if ((out=check_option(str,n,'m','u')) != NULL)
    sscanf(out,"%u",&mindim);
  if ((out=check_option(str,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(str,n,'r','f')) != NULL) {
    eps0set=1;
    sscanf(out,"%lf",&epsmin);
  }
  if ((out=check_option(str,n,'R','f')) != NULL) {
    eps1set=1;
    sscanf(out,"%lf",&epsmax);
  }
  if ((out=check_option(str,n,'#','u')) != NULL)
    sscanf(out,"%u",&epscount);
  if ((out=check_option(str,n,'n','u')) != NULL)
    sscanf(out,"%lu",&reference);
  if ((out=check_option(str,n,'s','u')) != NULL)
    sscanf(out,"%u",&maxiter);
  if ((out=check_option(str,n,'t','u')) != NULL)
    sscanf(out,"%u",&window);
  if ((out=check_option(str,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(str,n,'o','o')) != NULL)
    if (strlen(out) > 0)
      outfile=out;
}

void put_in_boxes(double eps)
{
  unsigned long i;
  long j,k;
  static unsigned long blength;

  blength=length-(maxdim-1)*delay-maxiter;

  for (i=0;i<BOX;i++)
    for (j=0;j<BOX;j++)
      box[i][j]= -1;

  for (i=0;i<blength;i++) {
    j=(long)(series[i]/eps)&ibox;
    k=(long)(series[i+delay]/eps)&ibox;
    liste[i]=box[j][k];
    box[j][k]=i;
  }
}

void lfind_neighbors(long act,double eps)
{
  unsigned int hi,k,k1;
  long i,j,i1,i2,j1,element;
  static long lwindow;
  double dx,eps2=sqr(eps);

  lwindow=(long)window;
  for (hi=0;hi<maxdim-1;hi++)
    found[hi]=0;
  i=(long)(series[act]/eps)&ibox;
  j=(long)(series[act+delay]/eps)&ibox;
  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1&ibox;
    for (j1=j-1;j1<=j+1;j1++) {
      element=box[i2][j1&ibox];
      while (element != -1) {
	if ((element < (act-lwindow)) || (element > (act+lwindow))) {
	  dx=sqr(series[act]-series[element]);
	  if (dx <= eps2) {
	    for (k=1;k<maxdim;k++) {
	      k1=k*delay;
	      dx += sqr(series[act+k1]-series[element+k1]);
	      if (dx <= eps2) {
		k1=k-1;
		lfound[k1][found[k1]]=element;
		found[k1]++;
	      }
	      else
		break;
	    }
	  }
	}
	element=liste[element];
      }
    }
  }
}

void iterate_points(long act)
{
  double **lfactor;
  double *dx;
  unsigned int i,j,l,l1;
  long k,element,**lcount;
  
  check_alloc(lfactor=(double**)malloc(sizeof(double*)*(maxdim-1)));
  check_alloc(lcount=(long**)malloc(sizeof(long*)*(maxdim-1)));
  for (i=0;i<maxdim-1;i++) {
    check_alloc(lfactor[i]=(double*)malloc(sizeof(double)*(maxiter+1)));
    check_alloc(lcount[i]=(long*)malloc(sizeof(long)*(maxiter+1)));
  }
  check_alloc(dx=(double*)malloc(sizeof(double)*(maxiter+1)));

  for (i=0;i<=maxiter;i++)
    for (j=0;j<maxdim-1;j++) {
      lfactor[j][i]=0.0;
      lcount[j][i]=0;
    }
  
  for (j=mindim-2;j<maxdim-1;j++) {
    for (k=0;k<found[j];k++) {
      element=lfound[j][k];
      for (i=0;i<=maxiter;i++)
	dx[i]=sqr(series[act+i]-series[element+i]);
      for (l=1;l<j+2;l++) {
	l1=l*delay;
	for (i=0;i<=maxiter;i++)
	  dx[i] += sqr(series[act+i+l1]-series[element+l1+i]);
      }
      for (i=0;i<=maxiter;i++)
	if (dx[i] > 0.0){
	  lcount[j][i]++;
	  lfactor[j][i] += dx[i];
	}
    }
  }
  for (i=mindim-2;i<maxdim-1;i++)
    for (j=0;j<=maxiter;j++)
      if (lcount[i][j]) {
	count[i][j]++;
	lyap[i][j] += log(lfactor[i][j]/lcount[i][j])/2.0;
      }
  
  for (i=0;i<maxdim-1;i++){
    free(lfactor[i]);
    free(lcount[i]);
  }
  free(lcount);
  free(lfactor);
  free(dx);
}

int main(int argc,char **argv)
{
  char stdi=0;
  double eps_fak;
  double epsilon;
  unsigned int i,j,l;
  FILE *fout;

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
      check_alloc(outfile=(char*)calloc(strlen(infile)+6,1));
      sprintf(outfile,"%s.lyap",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc(11,1));
      sprintf(outfile,"stdin.lyap");
    }
  }
  test_outfile(outfile);

  series=get_series(infile,&length,exclude,column,verbosity);
  rescale_data(series,length,&min,&max);

  if (eps0set)
    epsmin /= max;
  if (eps1set)
    epsmax /= max;

  if (epsmin >= epsmax) {
    epsmax=epsmin;
    epscount=1;
  }
  
  if (reference > (length-maxiter-(maxdim-1)*delay))
    reference=length-maxiter-(maxdim-1)*delay;
  if ((maxiter+(maxdim-1)*delay) >= length) {
    fprintf(stderr,"Too few points to handle these parameters!\n");
    exit(LYAP_K__MAXITER_TOO_LARGE);
  }

  if (maxdim < 2)
    maxdim=2;
  if (mindim < 2)
    mindim=2;
  if (mindim > maxdim)
    maxdim=mindim;
  
  check_alloc(liste=(long*)malloc(sizeof(long)*(length)));
  check_alloc(found=(long*)malloc(sizeof(long)*(maxdim-1)));
  check_alloc(lfound=(long**)malloc(sizeof(long*)*(maxdim-1)));
  for (i=0;i<maxdim-1;i++)
    check_alloc(lfound[i]=(long*)malloc(sizeof(long)*(length)));
  check_alloc(count=(long**)malloc(sizeof(long*)*(maxdim-1)));
  for (i=0;i<maxdim-1;i++)
    check_alloc(count[i]=(long*)malloc(sizeof(long)*(maxiter+1)));
  check_alloc(lyap=(double**)malloc(sizeof(double*)*(maxdim-1)));
  for (i=0;i<maxdim-1;i++)
    check_alloc(lyap[i]=(double*)malloc(sizeof(double)*(maxiter+1)));

  if (epscount == 1)
    eps_fak=1.0;
  else
    eps_fak=pow(epsmax/epsmin,1.0/(double)(epscount-1));

  fout=fopen(outfile,"w");
  if (verbosity&VER_INPUT)
    fprintf(stderr,"Opened %s for writing\n",outfile);
  for (l=0;l<epscount;l++) {
    epsilon=epsmin*pow(eps_fak,(double)l);
    for (i=0;i<maxdim-1;i++)
      for (j=0;j<=maxiter;j++) {
	count[i][j]=0;
	lyap[i][j]=0.0;
      }
    put_in_boxes(epsilon);
    for (i=0;i<reference;i++) {
      lfind_neighbors(i,epsilon);
      iterate_points(i);
    }
    if (verbosity&VER_USR1)
      fprintf(stderr,"epsilon= %e\n",epsilon*max);
    for (i=mindim-2;i<maxdim-1;i++) {
      fprintf(fout,"#epsilon= %e  dim= %d\n",epsilon*max,i+2);
      for (j=0;j<=maxiter;j++)
	if (count[i][j])
	  fprintf(fout,"%d %e %ld\n",j,lyap[i][j]/count[i][j],count[i][j]);
      fprintf(fout,"\n");
    }
    fflush(fout);
  }
  fclose(fout);
  return 0;
}
