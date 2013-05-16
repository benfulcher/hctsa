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
/*Author: Rainer Hegger, last modified: Apr 25, 2002 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the maximal Lyapunov exponent; Rosenstein et al."

#define NMAX 256

char *outfile=NULL;
char *infile=NULL;
char epsset=0;
double *series,*lyap;
long box[NMAX][NMAX],*list;
unsigned int dim=2,delay=1,steps=10,mindist=0;
unsigned int column=1;
unsigned int verbosity=0xff;
const unsigned int nmax=NMAX-1;
unsigned long length=ULONG_MAX,exclude=0;
long *found;
double eps0=1.e-3,eps,epsinv;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of datapoints [default is whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default is 0]\n");
  fprintf(stderr,"\t-c column to read[default 1]\n");
  fprintf(stderr,"\t-m embedding dimension [default 2]\n");
  fprintf(stderr,"\t-d delay  [default 1]\n");
  fprintf(stderr,"\t-t time window to omit [default 0]\n");
  fprintf(stderr,"\t-r epsilon size to start with [default "
	  "(data interval)/1000]\n");
  fprintf(stderr,"\t-s # of iterations [default 10]\n");
  fprintf(stderr,"\t-o name of output file [default 'datafile'.ros]\n");
  fprintf(stderr,"\t-V verbosity level [default 3]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
          "2='+ give more detailed information about the length scales\n");
  fprintf(stderr,"\t-h show these options\n");
  fprintf(stderr,"\n");
  exit(0);
}

void scan_options(int n,char **argv)
{
  char *out;

  if ((out=check_option(argv,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(argv,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(argv,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(argv,n,'m','u')) != NULL)
    sscanf(out,"%u",&dim);
  if ((out=check_option(argv,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(argv,n,'t','u')) != NULL)
    sscanf(out,"%u",&mindist);
  if ((out=check_option(argv,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&eps0);
  }
  if ((out=check_option(argv,n,'s','u')) != NULL)
    sscanf(out,"%u",&steps);
  if ((out=check_option(argv,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,n,'o','o')) != NULL)
    if (strlen(out) > 0)
      outfile=out;
}
      
void put_in_boxes(void)
{
  int i,j,x,y,del;
  
  for (i=0;i<NMAX;i++)
    for (j=0;j<NMAX;j++)
      box[i][j]= -1;

  del=delay*(dim-1);
  for (i=0;i<length-del-steps;i++) {
    x=(int)(series[i]*epsinv)&nmax;
    y=(int)(series[i+del]*epsinv)&nmax;
    list[i]=box[x][y];
    box[x][y]=i;
  }
}

char make_iterate(long act)
{
  char ok=0;
  int x,y,i,j,i1,k,del1=dim*delay;
  long element,minelement= -1;
  double dx,mindx=1.0;

  x=(int)(series[act]*epsinv)&nmax;
  y=(int)(series[act+delay*(dim-1)]*epsinv)&nmax;
  for (i=x-1;i<=x+1;i++) {
    i1=i&nmax;
    for (j=y-1;j<=y+1;j++) {
      element=box[i1][j&nmax];
      while (element != -1) {
	if (labs(act-element) > mindist) {
	  dx=0.0;
	  for (k=0;k<del1;k+=delay) {
	    dx += (series[act+k]-series[element+k])*
	      (series[act+k]-series[element+k]);
	    if (dx > eps*eps)
	      break;
	  }
	  if (k==del1) {
	    if (dx < mindx) {
	      ok=1;
	      if (dx > 0.0) {
		mindx=dx;
		minelement=element;
	      }
	    }
	  }
	}
	element=list[element];
      }
    }
  }
  if ((minelement != -1) ) {
    act--;
    minelement--;
    for (i=0;i<=steps;i++) {
      act++;
      minelement++;
      dx=0.0;
      for (j=0;j<del1;j+=delay) {
	dx += (series[act+j]-series[minelement+j])*
	  (series[act+j]-series[minelement+j]);
      }
      if (dx > 0.0) {
	found[i]++;
	lyap[i] += log(dx);
      }
    }
  }
  return ok;
}

int main(int argc,char **argv)
{
  char stdi=0,*done,alldone;
  int i;
  long n;
  long maxlength;
  double min,max;
  FILE *file;
  
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
      strcat(outfile,".ros");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.ros");
    }
  }
  test_outfile(outfile);

  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  rescale_data(series,length,&min,&max);

  if (epsset)
    eps0 /= max;

  check_alloc(list=(long*)malloc(length*sizeof(long)));
  check_alloc(lyap=(double*)malloc((steps+1)*sizeof(double)));
  check_alloc(found=(long*)malloc((steps+1)*sizeof(long)));
  check_alloc(done=(char*)malloc(length));

  for (i=0;i<=steps;i++) {
    lyap[i]=0.0;
    found[i]=0;
  }
  for (i=0;i<length;i++)
    done[i]=0;
  
  maxlength=length-delay*(dim-1)-steps-1-mindist;
  alldone=0;
  file=fopen(outfile,"w");
  if (verbosity&VER_INPUT)
    fprintf(stderr,"Opened %s for writing\n",outfile);
  for (eps=eps0;!alldone;eps*=1.1) {
    epsinv=1.0/eps;
    put_in_boxes();
    alldone=1;
    for (n=0;n<=maxlength;n++) {
      if (!done[n])
	done[n]=make_iterate(n);
      alldone &= done[n];
    }
    if (verbosity&VER_USR1)
      fprintf(stderr,"epsilon: %e already found: %ld\n",eps*max,found[0]);
  } 
  for (i=0;i<=steps;i++)
    if (found[i])
      fprintf(file,"%d %e\n",i,lyap[i]/found[i]/2.0);
  fclose(file);

  return 0;
}
