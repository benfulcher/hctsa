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
/*Author: Rainer Hegger. Last modified: Sep 29, 2000 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "routines/tsa.h"

#define WID_STR "Makes a local linear fit for multivariate data\n\
and iterates a trajectory"

#define NMAX 128

char onscreen=1,epsset=0,*outfile=NULL;
char *infile=NULL;
unsigned int nmax=(NMAX-1);
unsigned int verbosity=0xff;
long **box,*list,*found;
double **series,**cast;
double *interval,*min,epsilon;

unsigned int embed=2,dim=1,dim1,DELAY=1;
char *column=NULL,dimset=0,do_zeroth=0;
int MINN=30;
unsigned long LENGTH=ULONG_MAX,FLENGTH=1000,exclude=0;
double EPS0=1.e-3,EPSF=1.2;

double **mat,**imat,*vec,*localav,*foreav;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [Options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to be used [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default 0]\n");
  fprintf(stderr,"\t-c column [default 1,...,# of components]\n");
  fprintf(stderr,"\t-m #of components,embedding dimension [default 1,2]\n");
  fprintf(stderr,"\t-d delay for the embedding [default 1]\n");
  fprintf(stderr,"\t-L # of iterations [default 1000]\n");
  fprintf(stderr,"\t-k # of neighbors  [default 30]\n");
  fprintf(stderr,"\t-r size of initial neighborhood ["
	  " default (data interval)/1000]\n");
  fprintf(stderr,"\t-f factor to increase size [default 1.2]\n");
  fprintf(stderr,"\t-0 perfom a zeroth order fit [default not set]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.cast;"
	  " no -o means write to stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h  show these options\n");
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
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(in,n,'L','u')) != NULL)
    sscanf(out,"%lu",&FLENGTH);
  if ((out=check_option(in,n,'k','u')) != NULL)
    sscanf(out,"%u",&MINN);
  if ((out=check_option(in,n,'0','n')) != NULL)
    do_zeroth=1;
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&EPS0);
  }
  if ((out=check_option(in,n,'f','f')) != NULL)
    sscanf(out,"%lf",&EPSF);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    onscreen=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void put_in_boxes(void)
{
  int i,j,n;
  static int hdim;
  double epsinv;

  hdim=(embed-1)*DELAY;
  epsinv=1.0/epsilon;
  for (i=0;i<NMAX;i++)
    for (j=0;j<NMAX;j++)
      box[i][j]= -1;

  for (n=hdim;n<LENGTH-1;n++) {
    i=(int)(series[0][n]*epsinv)&nmax;
    j=(int)(series[dim1][n-hdim]*epsinv)&nmax;
    list[n]=box[i][j];
    box[i][j]=n;
  }
}

unsigned int hfind_neighbors(void)
{
  char toolarge;
  int i,j,i1,i2,j1,k,l,element;
  static int hdim;
  unsigned nfound=0;
  double max,dx,epsinv;

  hdim=(embed-1)*DELAY;
  epsinv=1.0/epsilon;
  i=(int)(cast[hdim][0]*epsinv)&nmax;
  j=(int)(cast[0][dim1]*epsinv)&nmax;
  
  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1&nmax;
    for (j1=j-1;j1<=j+1;j1++) {
      element=box[i2][j1&nmax];
      while (element != -1) {
	max=0.0;
	toolarge=0;
	for (l=0;l<dim;l++) {
	  for (k=0;k<=hdim;k += DELAY) {
	    dx=fabs(series[l][element-k]-cast[hdim-k][l]);
	    max=(dx>max) ? dx : max;
	    if (max > epsilon) {
	      toolarge=1;
	      break;
	    }
	  }
	  if (toolarge)
	    break;
	}
	if (max <= epsilon)
	  found[nfound++]=element;
	element=list[element];
      }
    }
  }
  return nfound;
}

void multiply_matrix(double **mat,double *vec)
{
  double *hvec;
  long i,j;

  check_alloc(hvec=(double*)malloc(sizeof(double)*dim*embed));
  for (i=0;i<dim*embed;i++) {
    hvec[i]=0.0;
    for (j=0;j<dim*embed;j++)
      hvec[i] += mat[i][j]*vec[j];
  }
  for (i=0;i<dim*embed;i++)
    vec[i]=hvec[i];
  free(hvec);
}

void make_fit(int number,double *newcast)
{
  double *sj,*si,lavi,lavj,fav;
  long i,i1,j,j1,hi,hj,hi1,hj1,n,which;
  static int hdim;

  hdim=(embed-1)*DELAY;

  for (i=0;i<dim*embed;i++)
    localav[i]=0.0;
  for (i=0;i<dim;i++)
    foreav[i]=0.0;

  for (n=0;n<number;n++) {
    which=found[n];
    for (j=0;j<dim;j++) {
      sj=series[j];
      foreav[j] += sj[which+1];
      for (j1=0;j1<embed;j1++) {
	hj=j*embed+j1;
	localav[hj] += sj[which-j1*DELAY];
      }
    }
  }

  for (i=0;i<dim*embed;i++)
    localav[i] /= number;
  for (i=0;i<dim;i++)
    foreav[i] /= number;

  for (i=0;i<dim;i++) {
    si=series[i];
    for (i1=0;i1<embed;i1++) {
      hi=i*embed+i1;
      lavi=localav[hi];
      hi1=i1*DELAY;
      for (j=0;j<dim;j++) {
	sj=series[j];
	for (j1=0;j1<embed;j1++) {
	  hj=j*embed+j1;
	  lavj=localav[hj];
	  hj1=j1*DELAY;
	  mat[hi][hj]=0.0;
	  if (hj >= hi) {
	    for (n=0;n<number;n++) {
	      which=found[n];
	      mat[hi][hj] += (si[which-hi1]-lavi)*(sj[which-hj1]-lavj);
	    }
	  }
	}
      }
    }
  }
  
  for (i=0;i<dim*embed;i++)
    for (j=i;j<dim*embed;j++) {
      mat[i][j] /= number;
      mat[j][i]=mat[i][j];
    }
  
  imat=invert_matrix(mat,dim*embed);

  for (i=0;i<dim;i++) {
    si=series[i];
    fav=foreav[i];
    for (j=0;j<dim;j++) {
      sj=series[j];
      for (j1=0;j1<embed;j1++) {
	hj=j*embed+j1;
	lavj=localav[hj];
	hj1=j1*DELAY;
	vec[hj]=0.0;
	for (n=0;n<number;n++) {
	  which=found[n];
	  vec[hj] += (si[which+1]-fav)*(sj[which-hj1]-lavj);
	}
	vec[hj] /= number;
      }
    }

    multiply_matrix(imat,vec);

    newcast[i]=foreav[i];
    for (j=0;j<dim;j++) {
      for (j1=0;j1<embed;j1++) {
	hj=j*embed+j1;
	newcast[i] += vec[hj]*(cast[hdim-j1*DELAY][j]-localav[hj]);
      }
    }
  }
  
  for (i=0;i<dim*embed;i++)
    free(imat[i]);
  free(imat);
}

void make_zeroth(int number,double *newcast)
{
  unsigned long i,d;
  double *sj;
  
  for (d=0;d<dim;d++) {
    newcast[d]=0.0;
    sj=series[d]+1;
    for (i=0;i<number;i++)
      newcast[d] += sj[found[i]];
    newcast[d] /= number;
  }
}

int main(int argc,char **argv)
{
  char stdi=0,done;
  long i,j,hdim,actfound;
  double maxinterval,*swap,*newcast;
  FILE *file=NULL;
  
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
      check_alloc(outfile=(char*)calloc(strlen(infile)+6,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".cast");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)11,(size_t)1));
      strcpy(outfile,"stdin.cast");
    }
  }
  if (!onscreen)
    test_outfile(outfile);
  
  hdim=(embed-1)*DELAY+1;
  if (column == NULL)
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dim,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&dim,column,
				      dimset,verbosity);
  check_alloc(min=(double*)malloc(sizeof(double)*dim));
  check_alloc(interval=(double*)malloc(sizeof(double)*dim));
  dim1=dim-1;
  maxinterval=0.0;
  for (i=0;i<dim;i++) {
    rescale_data(series[i],LENGTH,&min[i],&interval[i]);
    if (interval[i] > maxinterval)
      maxinterval=interval[i];
  }
  
  check_alloc(cast=(double**)malloc(sizeof(double*)*hdim));
  for (i=0;i<hdim;i++)
    check_alloc(cast[i]=(double*)malloc(sizeof(double)*dim));
  check_alloc(newcast=(double*)malloc(sizeof(double)*dim));
    
  check_alloc(list=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(found=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(box=(long**)malloc(sizeof(long*)*NMAX));
  for (i=0;i<NMAX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*NMAX));
  
  check_alloc(localav=(double*)malloc(sizeof(double)*dim*embed));
  check_alloc(foreav=(double*)malloc(sizeof(double)*dim));
  check_alloc(vec=(double*)malloc(sizeof(double)*dim*embed));
  check_alloc(mat=(double**)malloc(sizeof(double*)*dim*embed));
  for (i=0;i<dim*embed;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*dim*embed));

  if (epsset)
    EPS0 /= maxinterval;

  for (j=0;j<dim;j++)
    for (i=0;i<hdim;i++)
      cast[i][j]=series[j][LENGTH-hdim+i];
  
  if (!onscreen) {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  for (i=0;i<FLENGTH;i++) {
    done=0;
    epsilon=EPS0/EPSF;
    while (!done) {
      epsilon*=EPSF;
      put_in_boxes();
      actfound=hfind_neighbors();
      if (actfound >= MINN) {
	if (!do_zeroth)
	  make_fit(actfound,newcast);
	else
	  make_zeroth(actfound,newcast);
	if (onscreen) {
	  for (j=0;j<dim-1;j++)
	    printf("%e ",newcast[j]*interval[j]+min[j]);
	  printf("%e\n",newcast[dim-1]*interval[dim-1]+min[dim-1]);
	  fflush(stdout);
	}
	else {
	  for (j=0;j<dim-1;j++)
	    fprintf(file,"%e ",newcast[j]*interval[j]+min[j]);
	  fprintf(file,"%e\n",newcast[dim-1]*interval[dim-1]+min[dim-1]);
	  fflush(file);
	}
	done=1;
	for (j=0;j<dim;j++) {
	  if ((newcast[j] > 2.0) || (newcast[j] < -1.0)) {
	    fprintf(stderr,"Forecast failed. Escaping data region!\n");
	    exit(NSTEP__ESCAPE_REGION);
	  }
	}

	swap=cast[0];
	for (j=0;j<hdim-1;j++)
	  cast[j]=cast[j+1];
	cast[hdim-1]=swap;
	for (j=0;j<dim;j++)
	  cast[hdim-1][j]=newcast[j];
      }
    }
  }
  if (!onscreen)
    fclose(file);
  
  if (outfile != NULL)
    free(outfile);
  for (i=0;i<embed*dim;i++)
    free(mat[i]);
  free(mat);
  for (i=0;i<hdim;i++)
    free(cast[i]);
  free(cast);
  free(newcast);
  free(found);
  free(list);
  for (i=0;i<NMAX;i++)
    free(box[i]);
  free(box);
  free(vec);
  free(localav);
  free(foreav);
  free(min);
  free(interval);
  for (i=0;i<dim;i++)
    free(series[i]);
  free(series);

  return 0;
}
