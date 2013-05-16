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
/*Author: Rainer Hegger, last modified: Mar 20, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Estimates the finite size Lyapunov exponent; Vulpiani et al."


#define NMAX 256

char *outfile=NULL;
char *infile=NULL;
char epsset=0,stdo=1;
double *series;
long box[NMAX][NMAX],*list;
unsigned int dim=2,delay=1,mindist=0;
unsigned int column=1;
unsigned int verbosity=0xff;
const unsigned int nmax=NMAX-1;
unsigned long length=ULONG_MAX,exclude=0;
double eps0=1.e-3,eps,epsinv,epsmax,epsfactor;
int howmany;

struct fsle {
  double time,factor,eps;
  long count;
} *data;

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
	  "(std. dev. of data)/1000]\n");
  fprintf(stderr,"\t-o name of output file [default 'datafile'.fsl ,"
	  "without -o: stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
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
  if ((out=check_option(argv,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&eps0);
  }
  if ((out=check_option(argv,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}
      
void put_in_boxes(void)
{
  int i,j,x,y,del;

  for (i=0;i<NMAX;i++)
    for (j=0;j<NMAX;j++)
      box[i][j]= -1;

  del=delay*(dim-1);
  for (i=0;i<length-del;i++) {
    x=(int)(series[i]*epsinv)&nmax;
    y=(int)(series[i+del]*epsinv)&nmax;
    list[i]=box[x][y];
    box[x][y]=i;
  }
}

char make_iterate(long act)
{
  char ok=0;
  int x,y,i,j,i1,k,del1=dim*delay,which;
  long element,minelement= -1;
  double dx=0.0,mindx=2.0,stime;

  x=(int)(series[act]*epsinv)&nmax;
  y=(int)(series[act+delay*(dim-1)]*epsinv)&nmax;
  for (i=x-1;i<=x+1;i++) {
    i1=i&nmax;
    for (j=y-1;j<=y+1;j++) {
      element=box[i1][j&nmax];
      while (element != -1) {
	if (labs(act-element) > mindist) {
	  for (k=0;k<del1;k+=delay) {
	    dx = fabs(series[act+k]-series[element+k]);
	    if (dx > eps)
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
  
  if ((minelement != -1) && (mindx < eps)) {
    act += del1-delay+1;
    minelement += del1-delay+1;
    which=(int)(log(mindx/eps0)/log(epsfactor));
    if (which < 0) {
      while ((dx=fabs(series[act]-series[minelement])) < data[0].eps) {
	act++;
	minelement++;
	if ((act >= length) || (minelement >= length))
	  return ok;
      }
      mindx=dx;
      which=(int)(log(mindx/eps0)/log(epsfactor));
    }
    for (i=which;i<howmany-1;i++) {
      stime=0;
      while ((dx=fabs(series[act]-series[minelement])) < data[i+1].eps) {
	act++;
	minelement++;
	if ((act >= length) || (minelement >= length))
	  return ok;
	stime++;
      }
      if (stime > 0) {
	data[i].time += stime;
	data[i].factor += log(dx/mindx);
	data[i].count++;
      }
      mindx=dx;
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
  double min,max,se_av,se_var,se0_av,se0_var;
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
      strcat(outfile,".fsl");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.fsl");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  variance(series,length,&se0_av,&se0_var);
  rescale_data(series,length,&min,&max);
  variance(series,length,&se_av,&se_var);
  
  if (epsset) {
    eps0 /= max;
    epsmax=se0_var;
  }
  else {
    eps0 *= se_var;
    epsmax=se_var;
  }
  if (eps0 >= epsmax) {
    fprintf(stderr,"The minimal epsilon is too large. Exiting!\n");
    exit(FSLE__TOO_LARGE_MINEPS);
  }
  epsfactor=sqrt(2.0);

  howmany=(int)(log(epsmax/eps0)/log(epsfactor))+1;
  check_alloc(data=(struct fsle*)malloc(sizeof(struct fsle)*howmany));
  eps=eps0/epsfactor;
  for (i=0;i<howmany;i++) {
    data[i].time=data[i].factor=0.0;
    data[i].eps= (eps *= epsfactor);
    data[i].count=0;
  }
  
  check_alloc(list=(long*)malloc(length*sizeof(long)));
  check_alloc(done=(char*)malloc(length));

  for (i=0;i<length;i++)
    done[i]=0;
  
  maxlength=length-delay*(dim-1)-1-mindist;
  alldone=0;
  for (eps=eps0;(eps<=epsmax) && (!alldone);eps*=epsfactor) {
    epsinv=1.0/eps;
    put_in_boxes();
    alldone=1;
    for (n=0;n<=maxlength;n++) {
      if (!done[n])
	done[n]=make_iterate(n);
      alldone &= done[n];
    }
  } 
  if (!stdo) {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    for (i=0;i<howmany;i++)
      if (data[i].factor > 0.0)
	fprintf(file,"%e %e %ld\n",data[i].eps*max,
		data[i].factor/data[i].time,data[i].count);
    fclose(file);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    for (i=0;i<howmany;i++)
      if (data[i].factor > 0.0)
	fprintf(stdout,"%e %e %ld\n",data[i].eps*max,
		data[i].factor/data[i].time,data[i].count);
  }    

  if (infile != NULL)
    free(infile);
  if (outfile != NULL)
    free(outfile);
  free(series);
  free(data);
  free(list);
  free(done);

  return 0;
}
