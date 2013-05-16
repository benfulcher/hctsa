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
/*Author: Rainer Hegger. Last modified: Feb 19, 2007 */
/* Changes:
     2/19/2007:  Changed name and default for noise  
     10/26/2006: Add seed option
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "routines/tsa.h"

#define WID_STR "Makes a local zeroth order forecast for multivariate data\n\
and iterates a trajectory"

#define NMAX 128

char onscreen=1,epsset=0,*outfile=NULL,setsort=1,setnoise=0;
char *infile=NULL;
unsigned int nmax=(NMAX-1);
unsigned int verbosity=0xff;
long **box,*list,*found;
double **series,**cast,*abstand,*var;
double epsilon;

unsigned int embed=2,dim=1,dim1,DELAY=1;
char *column=NULL,dimset=0;
unsigned int MINN=50;
unsigned int **indexes;
unsigned long LENGTH=ULONG_MAX,FLENGTH=1000,exclude=0;
unsigned long seed=0x9074325L;
double EPS0=1.e-3,EPSF=1.2,Q=10.0;

double **mat,*vec,*hsum,*newav;

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
  fprintf(stderr,"\t-k # of neighbors  [default %u]\n",MINN);
  fprintf(stderr,"\t-K fix # of neighbors  [default no]\n");
  fprintf(stderr,"\t-%% # variance of noise [default %3.1lf]\n",Q);
  fprintf(stderr,"\t-I seed for the rnd-generator (If seed=0, the time\n"
          "\t\tcommand is used to set the seed) [Default: fixed]\n");
  fprintf(stderr,"\t-r size of initial neighborhood ["
	  " default (data interval)/1000]\n");
  fprintf(stderr,"\t-f factor to increase size [default 1.2]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.lzr;"
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
  if ((out=check_option(in,n,'K','n')) != NULL)
    setsort=1;
  if ((out=check_option(in,n,'I','u')) != NULL) {
    sscanf(out,"%lu",&seed);
    if (seed == 0)
      seed=(unsigned long)time((time_t*)&seed);
  }
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&EPS0);
  }
  if ((out=check_option(in,n,'f','f')) != NULL)
    sscanf(out,"%lf",&EPSF);
  if ((out=check_option(in,n,'%','f')) != NULL) {
    sscanf(out,"%lf",&Q);
    if (Q>0.0)
      setnoise=1;
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    onscreen=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void sort(unsigned long nfound)
{
  double dx,dswap;
  int i,j,k,hf,iswap,hdim;

  hdim=(embed-1)*DELAY;

  for (i=0;i<nfound;i++) {
    hf=found[i];
    abstand[i]=0.0;
    for (j=0;j<dim;j++) {
      for (k=0;k<=hdim;k += DELAY) {
	dx=fabs(series[j][hf-k]-cast[hdim-k][j]);
	if (dx > abstand[i]) abstand[i]=dx;
      }
    }
  }

  for (i=0;i<MINN;i++)
    for (j=i+1;j<nfound;j++)
      if (abstand[j]<abstand[i]) {
	dswap=abstand[i];
	abstand[i]=abstand[j];
	abstand[j]=dswap;
	iswap=found[i];
	found[i]=found[j];
	found[j]=iswap;
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
  int i,j,i1,i2,j1,l,hc,hd,element;
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
	for (l=0;l<dim*embed;l++) {
	  hc=indexes[0][l];
	  hd=indexes[1][l];
	  dx=fabs(series[hc][element-hd]-cast[hdim-hd][hc]);
	  max=(dx>max) ? dx : max;
	  if (max > epsilon) {
	    toolarge=1;
	    break;
	  }
	}
	if (max <= epsilon)
	  found[nfound++]=element;
	element=list[element];
      }
    }
  }
  return nfound;
}

void make_zeroth(int number,double *newcast)
{
  long d,i;
  double *sd;
  
  for (d=0;d<dim;d++) {
    newcast[d]=0.0;
    sd=series[d]+1;
    for (i=0;i<number;i++)
      newcast[d] += sd[found[i]];
    newcast[d] /= (double)number;
  }

  if (setnoise) {
    for (d=0;d<dim;d++)
      newcast[d] += gaussian(var[d]*Q);
  }
}

int main(int argc,char **argv)
{
  char stdi=0,done;
  long i,j,hdim,actfound;
  unsigned long count=1;
  double *swap,*newcast,maxinterval,*min,*interval,dummy,epsilon0;
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
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".lzr");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.lzr");
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

  dim1=dim-1;

  check_alloc(min=(double*)malloc(sizeof(double)*dim));
  check_alloc(interval=(double*)malloc(sizeof(double)*dim));
  check_alloc(var=(double*)malloc(sizeof(double)*dim));

  maxinterval=0.0;

  for (i=0;i<dim;i++) {
    rescale_data(series[i],LENGTH,&min[i],&interval[i]);
    variance(series[i],LENGTH,&dummy,&var[i]);
    if (interval[i] > maxinterval)
      maxinterval=interval[i];
  }

  if (epsset)
    EPS0 /= maxinterval;
    
  check_alloc(cast=(double**)malloc(sizeof(double*)*hdim));
  for (i=0;i<hdim;i++)
    check_alloc(cast[i]=(double*)malloc(sizeof(double)*dim));
  check_alloc(newcast=(double*)malloc(sizeof(double)*dim));
  check_alloc(newav=(double*)malloc(sizeof(double)*dim));
    
  check_alloc(list=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(found=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(abstand=(double*)malloc(sizeof(double)*LENGTH));
  check_alloc(box=(long**)malloc(sizeof(long*)*NMAX));
  for (i=0;i<NMAX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*NMAX));
  
  check_alloc(vec=(double*)malloc(sizeof(double)*dim));
  check_alloc(hsum=(double*)malloc(sizeof(double)*dim));
  check_alloc(mat=(double**)malloc(sizeof(double*)*dim));
  for (i=0;i<dim;i++) {
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*dim));
  }

  for (j=0;j<dim;j++)
    for (i=0;i<hdim;i++)
      cast[i][j]=series[j][LENGTH-hdim+i];

  indexes=make_multi_index(dim,embed,DELAY);
  
  if (!onscreen) {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  rnd_init(seed);

  epsilon0=EPS0/EPSF;

  if (setnoise) 
    Q /= 100.0;

  for (i=0;i<FLENGTH;i++) {
    done=0;
    if (setsort)
      epsilon= epsilon0/((double)count*EPSF);
    else
      epsilon=epsilon0;
    while (!done) {
      epsilon*=EPSF;
      put_in_boxes();
      actfound=hfind_neighbors();
      if (actfound >= MINN) {
	if (setsort) {
	  epsilon0 += epsilon;
	  count++;
	  sort(actfound);
	  actfound=MINN;
	}
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
  for (i=0;i<dim;i++)
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
  free(newav);
  for (i=0;i<dim;i++)
    free(series[i]);
  free(series);

  return 0;
}
