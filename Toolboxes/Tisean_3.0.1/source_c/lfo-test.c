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
/*Author: Rainer Hegger */
/*Changes:
  Sep 8, 2006: Add -o functionality
  Sep 7, 2006: Completely rewritten to handle multivariate data
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"
#include <math.h>

#define WID_STR "Estimates the average forecast error of a local\n\t\
linear fit"


/*number of boxes for the neighbor search algorithm*/
#define NMAX 512

unsigned int nmax=(NMAX-1),comp1,hdim,**indexes;
long **box,*list;
unsigned long *found,*hfound;
double **series;
double epsilon;
double **mat,**imat,*vec,*localav,*foreav;

char epsset=0,causalset=0;
unsigned int verbosity=VER_INPUT|VER_FIRST_LINE;
unsigned int COMP=1,EMBED=2,DIM,DELAY=1,MINN=30,STEP=1;
double EPS0=1.e-3,EPSF=1.2;
unsigned long LENGTH=ULONG_MAX,exclude=0,CLENGTH=ULONG_MAX,causal;
char *infile=NULL,*COLUMN=NULL,*outfile=NULL;
char dimset=0,stout=1;

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
  fprintf(stderr,"\t-c columns to read [default: 1]\n");
  fprintf(stderr,"\t-m # of components, embedding dimension "
	  "[default: %u,%u]\n",COMP,EMBED);
  fprintf(stderr,"\t-d delay [default: 1]\n");
  fprintf(stderr,"\t-n iterations [default: length]\n");
  fprintf(stderr,"\t-k minimal number of neighbors for the fit "
	  "[default: 30]\n");
  fprintf(stderr,"\t-r neighborhoud size to start with "
	  "[default: (data interval)/1000]\n");
  fprintf(stderr,"\t-f factor to increase size [default: 1.2]\n");
  fprintf(stderr,"\t-s steps to forecast [default: 1]\n");
  fprintf(stderr,"\t-C width of causality window [default: steps]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.fce"
	  " no -o means write to stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
	  "2='+ print indiviual forecast errors'\n");
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
    COLUMN=out;
    dimset=1;
  }
  if ((out=check_option(in,n,'m','2')) != NULL)
    sscanf(out,"%u,%u",&COMP,&EMBED);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(in,n,'n','u')) != NULL)
    sscanf(out,"%lu",&CLENGTH);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
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
  if ((out=check_option(in,n,'C','u')) != NULL) {
    sscanf(out,"%lu",&causal);
    causalset=1;
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void put_in_boxes(void)
{
  int i,j,n;
  double epsinv;

  epsinv=1.0/epsilon;
  for (i=0;i<NMAX;i++)
    for (j=0;j<NMAX;j++)
      box[i][j]= -1;

  for (n=hdim;n<LENGTH-STEP;n++) {
    i=(int)(series[0][n]*epsinv)&nmax;
    j=(int)(series[comp1][n-hdim]*epsinv)&nmax;
    list[n]=box[i][j];
    box[i][j]=n;
  }
}

unsigned int hfind_neighbors(unsigned long act)
{
  char toolarge;
  int i,j,i1,i2,j1,k,element;
  unsigned long nfound=0;
  unsigned int hcomp,hdel;
  double max,dx,epsinv;

  epsinv=1.0/epsilon;

  i=(int)(series[0][act]*epsinv)&nmax;
  j=(int)(series[comp1][act-hdim]*epsinv)&nmax;
  
  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1&nmax;
    for (j1=j-1;j1<=j+1;j1++) {
      element=box[i2][j1&nmax];
      while (element != -1) {
	max=0.0;
	toolarge=0;
	for (k=0;k<DIM;k += 1) {
	  hcomp=indexes[0][k];
	  hdel=indexes[1][k];
	  dx=fabs(series[hcomp][element-hdel]-series[hcomp][act-hdel]);
	  max=(dx>max) ? dx : max;
	  if (max > epsilon) {
	    toolarge=1;
	    break;
	  }
	  if (toolarge)
	    break;
	}
	if (max <= epsilon)
	  hfound[nfound++]=element;
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

  check_alloc(hvec=(double*)malloc(sizeof(double)*DIM));
  for (i=0;i<DIM;i++) {
    hvec[i]=0.0;
    for (j=0;j<DIM;j++)
      hvec[i] += mat[i][j]*vec[j];
  }
  for (i=0;i<DIM;i++)
    vec[i]=hvec[i];
  free(hvec);
}

void make_fit(int number,unsigned long act,double *newcast)
{
  double *sj,*si,lavi,lavj,fav;
  unsigned int hci,hdi,hcj,hdj;
  long i,j,n,which;

  for (i=0;i<DIM;i++)
    localav[i]=0.0;
  for (i=0;i<COMP;i++)
    foreav[i]=0.0;

  for (n=0;n<number;n++) {
    which=found[n];
    for (j=0;j<COMP;j++)
      foreav[j] += series[j][which+STEP];
    for (j=0;j<DIM;j++) {
      hcj=indexes[0][j];
      hdj=indexes[1][j];
      localav[j] += series[hcj][which-hdj];
    }
  }

  for (i=0;i<DIM;i++)
    localav[i] /= number;
  for (i=0;i<COMP;i++)
    foreav[i] /= number;

  for (i=0;i<DIM;i++) {
    hci=indexes[0][i];
    hdi=indexes[1][i];
    lavi=localav[i];
    si=series[hci];
    for (j=i;j<DIM;j++) {
      hcj=indexes[0][j];
      hdj=indexes[1][j];
      lavj=localav[j];
      sj=series[hcj];
      mat[i][j]=0.0;
      for (n=0;n<number;n++) {
	which=found[n];
	mat[i][j] += (si[which-hdi]-lavi)*(sj[which-hdj]-lavj);
      }
      mat[i][j] /= number;
      mat[j][i] = mat[i][j];
    }
  }

  imat=invert_matrix(mat,DIM);

  for (i=0;i<COMP;i++) {
    si=series[i];
    fav=foreav[i];
    for (j=0;j<DIM;j++) {
      hcj=indexes[0][j];
      hdj=indexes[1][j];
      lavj=localav[j];
      vec[j]=0.0;
      sj=series[hcj];
      for (n=0;n<number;n++) {
	which=found[n];
	vec[j] += (si[which+STEP]-fav)*(sj[which-hdj]);
      }
      vec[j] /= number;
    }

    multiply_matrix(imat,vec);

    newcast[i]=foreav[i];
    for (j=0;j<DIM;j++) {
      hcj=indexes[0][j];
      hdj=indexes[1][j];
      newcast[i] += vec[j]*(series[hcj][act-hdj]-localav[j]);
    }
  }
  

  for (i=0;i<DIM;i++)
    free(imat[i]);
  free(imat);
}

int main(int argc,char **argv)
{
  char stin=0,alldone,*done;
  long i,j;
  unsigned long actfound;
  unsigned long clength;
  double *rms,*av,*min,*interval,maxinterval,norm;
  double *error,**individual=NULL;
  double *newcast;
  FILE *fout;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);

  if (!causalset)
    causal=STEP;

#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stin=1;
  
  if (outfile == NULL) {
    if (!stin) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".fce");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.fce");
    }
  }
  if (!stout)
    test_outfile(outfile);
  
  if (COLUMN == NULL)
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&COMP,"",dimset,
                                      verbosity);
  else
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&COMP,COLUMN,
                                      dimset,verbosity);

  if ((LENGTH-(EMBED-1)*DELAY) < MINN) {
    fprintf(stderr,"Data set is too short to find enough neighbors "
	    "for the fit! Exiting!\n");
    exit(ONESTEP_TOO_FEW_POINTS);
  }

  DIM=EMBED*COMP;
  check_alloc(min=(double*)malloc(sizeof(double)*COMP));
  check_alloc(interval=(double*)malloc(sizeof(double)*COMP));
  check_alloc(av=(double*)malloc(sizeof(double)*COMP));
  check_alloc(rms=(double*)malloc(sizeof(double)*COMP));

  maxinterval=0.0;
  for (i=0;i<COMP;i++) {
    rescale_data(series[i],LENGTH,&min[i],&interval[i]);
    maxinterval=(maxinterval<interval[i])?interval[i]:maxinterval;
    variance(series[i],LENGTH,&av[i],&rms[i]);
  }
  
  if (verbosity&VER_USR1) {
    check_alloc(individual=(double**)malloc(sizeof(double*)*COMP));
    for (j=0;j<COMP;j++) {
      check_alloc(individual[j]=(double*)malloc(sizeof(double)*LENGTH));
      for (i=0;i<LENGTH;i++)
	individual[j][i]=0.0;
    }
  }

  check_alloc(list=(long*)malloc(sizeof(long)*LENGTH));
  check_alloc(found=(unsigned long*)malloc(sizeof(long)*LENGTH));
  check_alloc(hfound=(unsigned long*)malloc(sizeof(long)*LENGTH));
  check_alloc(done=(char*)malloc(sizeof(char)*LENGTH));
  check_alloc(box=(long**)malloc(sizeof(long*)*NMAX));
  for (i=0;i<NMAX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*NMAX));
    
  for (i=0;i<LENGTH;i++)
    done[i]=0;

  alldone=0;
  if (epsset)
    EPS0 /= maxinterval;

  epsilon=EPS0/EPSF;
  clength=(CLENGTH <= LENGTH) ? CLENGTH-STEP : LENGTH-STEP;
  comp1=COMP-1;
  indexes=make_multi_index(COMP,EMBED,DELAY);

  hdim=(EMBED-1)*DELAY;
  check_alloc(newcast=(double*)malloc(sizeof(double)*COMP));


  check_alloc(localav=(double*)malloc(sizeof(double)*DIM));
  check_alloc(foreav=(double*)malloc(sizeof(double)*COMP));
  check_alloc(vec=(double*)malloc(sizeof(double)*DIM));
  check_alloc(mat=(double**)malloc(sizeof(double*)*DIM));
  for (i=0;i<=DIM;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*DIM));

  check_alloc(error=(double*)malloc(sizeof(double)*COMP));
  for (i=0;i<COMP;i++)
    error[i]=0.0;

  while (!alldone) {
    alldone=1;
    epsilon*=EPSF;
    put_in_boxes() ;
    for (i=(EMBED-1)*DELAY;i<clength;i++)
      if (!done[i]) {
	actfound=hfind_neighbors(i);
	actfound=exclude_interval(actfound,i-causal+1,
				  i+causal+(EMBED-1)*DELAY-1,hfound,found);
	if (actfound > MINN) {
	  make_fit(actfound,i,newcast);
	  for (j=0;j<COMP;j++)
	    error[j] += sqr(newcast[j]-series[j][i+STEP]);
	  if (verbosity&VER_USR1) {
	    for (j=0;j<COMP;j++)
	      individual[j][i]=(newcast[j]-series[j][i+STEP])*interval[j];
	  }
	  done[i]=1;
	}
	alldone &= done[i];
      }
  }
  norm=((double)clength-(double)((EMBED-1)*DELAY));
  if (stout) {
    if (verbosity&VER_USR1) {
      fprintf(stdout,"#Relative forecast errors for each component:\n");
      for (i=0;i<COMP;i++) 
	fprintf(stdout,"# %e\n",sqrt(error[i]/norm)/rms[i]);
    
      for (i=(EMBED-1)*DELAY;i<clength;i++) {
	for (j=0;j<COMP-1;j++)
	  fprintf(stdout,"%e ",individual[j][i]);
	fprintf(stdout,"%e\n",individual[COMP-1][i]);
      }
    }
    else {
      fprintf(stdout,"#Relative forecast errors for each component:\n");
      for (i=0;i<COMP;i++) 
	fprintf(stdout,"%e\n",sqrt(error[i]/norm)/rms[i]);
    }
  }
  else {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    if (verbosity&VER_USR1) {
      fprintf(fout,"#Relative forecast errors for each component:\n");
      for (i=0;i<COMP;i++) 
	fprintf(fout,"# %e\n",sqrt(error[i]/norm)/rms[i]);
    
      for (i=(EMBED-1)*DELAY;i<clength;i++) {
	for (j=0;j<COMP-1;j++)
	  fprintf(fout,"%e ",individual[j][i]);
	fprintf(fout,"%e\n",individual[COMP-1][i]);
      }
    }
    else {
      fprintf(fout,"#Relative forecast errors for each component:\n");
      for (i=0;i<COMP;i++) 
	fprintf(fout,"%e\n",sqrt(error[i]/norm)/rms[i]);
    }
    fclose(fout);
    free(outfile);
  }

  return 0;
}
