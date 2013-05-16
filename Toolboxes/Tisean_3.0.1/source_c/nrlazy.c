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
/*Author: Rainer Hegger Last modified: Nov 30, 2000 */
/*Changes:
  12/11/05: Going multivariate
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Performs simple noise reduction."

#define BOX (unsigned int)512

unsigned long length=ULONG_MAX,exclude=0;
unsigned int comp=1,embed=5,delay=1,iterations=1,alldim;
unsigned int verbosity=0x3;
char *column=NULL;
double eps=1.0e-3,epsvar;

char *outfile=NULL,epsset=0,stdo=1,epsvarset=0;
char *infile=NULL;
double **series,**corr,*interval,*min,*hcor;
long **box,*list,**nf;
unsigned int **indexes;
char dimset=0;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [Options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to use [default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default: 0]\n");
  fprintf(stderr,"\t-c column to read [default: 1]\n");
  fprintf(stderr,"\t-m no. of comp.,embedding dim. [default: %u,%u]\n",
	  comp,embed);
  fprintf(stderr,"\t-d delay [default: 1]\n");
  fprintf(stderr,"\t-i iterations [default: 1]\n");
  fprintf(stderr,"\t-r neighborhoud size [default: (interval of data)/1000]\n");
  fprintf(stderr,"\t-v neighborhoud size (in units of the std. dev. of the "
	  "data \n\t\t(overwrites -r) [default: not set]\n");
  fprintf(stderr,"\t-o output file name [Default: 'datafile'.laz.n,"
	  "\n\t\twhere n is the number of the last iteration,"
	  "\n\t\twithout -o the last iteration is written to stdout.]\n");
  fprintf(stderr,"\t-V verbosity level [Default: 3]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
	  "2='+ write output of all iterations to files'\n\t\t"
	  "4='+ write the number of neighbors found for each point\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **in)
{
  char *out;

  if ((out=check_option(in,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'c','c')) != NULL) {
    column=out;
    dimset=1;
  }
  if ((out=check_option(in,n,'m','2')) != NULL)
    sscanf(out,"%u,%u",&comp,&embed);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(in,n,'i','u')) != NULL)
    sscanf(out,"%u",&iterations);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&eps);
  }
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'v','f')) != NULL) {
    epsvarset=1;
    sscanf(out,"%lf",&epsvar);
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

unsigned int correct(unsigned long n)
{
  int i,i1,i2,j,j1,k;
  int ibox=BOX-1;
  unsigned int hdel,hcomp;
  double epsinv,dx;
  long element,nfound=0;

  epsinv=1./eps;

  for (i=0;i<alldim;i++)
    hcor[i]=0.0;

  i=(int)(series[0][n]*epsinv)&ibox;
  j=(int)(series[comp-1][n-(embed-1)*delay]*epsinv)&ibox;
  
  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1&ibox;
    for (j1=j-1;j1<=j+1;j1++) {
      element=box[i2][j1&ibox];
      while (element != -1) {
	for (k=0;k<alldim;k++) {
	  hcomp=indexes[0][k];
	  hdel=indexes[1][k];
	  dx=fabs(series[hcomp][n-hdel]-series[hcomp][element-hdel]);
	  if (dx > eps)
	    break;
	}
	if (k == alldim) {
	  nfound++;
	  for (k=0;k<alldim;k++) {
	    hcomp=indexes[0][k];
	    hdel=indexes[1][k];
	    hcor[k] += series[hcomp][element-hdel];
	  }
	}
	element=list[element];
      }
    }
  }
  for (k=0;k<alldim;k++) {
    hcomp=indexes[0][k];
    hdel=indexes[1][k];
    corr[hcomp][n-hdel] += hcor[k]/nfound;
    nf[hcomp][n-hdel]++;
  }

  return nfound;
}

int main(int argc,char **argv)
{
  char *ofname;
  char stdi=0;
  int iter;
  unsigned int *nmf;
  unsigned long n,i;
  double dav,dvar,maxinterval,maxdvar;
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
      check_alloc(ofname=(char*)calloc(strlen(infile)+9,(size_t)1));
      sprintf(outfile,"%s.laz",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      check_alloc(ofname=(char*)calloc((size_t)14,(size_t)1));
      sprintf(outfile,"stdin.laz");
    }
  }
  else
    check_alloc(ofname=(char*)calloc(strlen(outfile)+10,(size_t)1));


  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&comp,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&comp,column,
				      dimset,verbosity);

  check_alloc(interval=(double*)malloc(sizeof(double)*comp));
  check_alloc(min=(double*)malloc(sizeof(double)*comp));

  maxinterval=maxdvar=0.0;
  for (i=0;i<comp;i++) {
    rescale_data(series[i],length,&min[i],&interval[i]);
    if (interval[i] > maxinterval) maxinterval=interval[i];
    variance(series[i],length,&dav,&dvar);
    if (dvar > maxdvar)  maxdvar=dvar;
  }
  alldim=comp*embed;

  check_alloc(nmf=(unsigned int*)malloc(sizeof(int)*length));
  check_alloc(list=(long*)malloc(sizeof(long)*length));
  check_alloc(box=(long**)malloc(sizeof(long*)*BOX));
  for (n=0;n<BOX;n++)
    check_alloc(box[n]=(long*)malloc(sizeof(long)*BOX));

  check_alloc(nf=(long**)malloc(sizeof(long*)*comp));
  check_alloc(corr=(double**)malloc(sizeof(double*)*comp));
  for (i=0;i<comp;i++) {
    check_alloc(nf[i]=(long*)malloc(sizeof(long)*length));
    check_alloc(corr[i]=(double*)malloc(sizeof(double)*length));
  }

  indexes=make_multi_index(comp,embed,delay);

  if (epsset)
    eps/=maxinterval;
  else
    eps=1.0/1000.;

  if (epsvarset)
    eps=epsvar*maxdvar;

  for (iter=1;iter<=iterations;iter++) {
    make_multi_box2(series,box,list,length,BOX,comp,embed,delay,eps);
    for (n=0;n<length;n++) {
      for (i=0;i<comp;i++) {
	corr[i][n]=0.0;
	nf[i][n]=0;
      }
      nmf[n]=1;
    }
    
    check_alloc(hcor=(double*)malloc(sizeof(double)*alldim));
    for (n=(embed-1)*delay;n<length;n++)
      nmf[n]=correct(n);
    free(hcor);
    
    for (n=0;n<length;n++)
      for (i=0;i<comp;i++)
	if (nf[i][n])
	  series[i][n]=corr[i][n]/nf[i][n];

    if ((verbosity&VER_USR1) && (iter < iterations)) {
      sprintf(ofname,"%s.%d",outfile,iter);
      test_outfile(ofname);
      file=fopen(ofname,"w");
      if (verbosity&VER_INPUT)
	fprintf(stderr,"Opened %s for writing\n",ofname);
      if (stdo && (iter == iterations)) {
	if (verbosity&VER_INPUT)
	  fprintf(stderr,"Writing to stdout\n");
      }
      for (n=0;n<length;n++) {
	if (stdo && (iter == iterations)) {
	  if (verbosity&VER_USR2) {
	    for (i=0;i<comp;i++) 
	      fprintf(stdout,"%e ",series[i][n]*interval[i]+min[i]);
	    fprintf(stdout,"%u\n",nmf[n]);
	  }
	  else {
	    fprintf(stdout,"%e",series[0][n]*interval[0]+min[0]);
	    for (i=1;i<comp;i++)
	      fprintf(stdout,"%e ",series[i][n]*interval[i]+min[i]);
	    fprintf(stdout,"\n");
	  }
	}
	if (verbosity&VER_USR2) {
	  for (i=0;i<comp;i++) 
	    fprintf(file,"%e ",series[i][n]*interval[i]+min[i]);
	  fprintf(file,"%u\n",nmf[n]);
	}
	else {
	  fprintf(file,"%e",series[0][n]*interval[0]+min[0]);
	  for (i=1;i<comp;i++)
	    fprintf(file," %e",series[i][n]*interval[i]+min[i]);
	  fprintf(file,"\n");
	}
      }
      fclose(file);
    }
    if (iter == iterations) {
      if (!stdo || (verbosity&VER_USR1)) {
	sprintf(ofname,"%s.%d",outfile,iter);
	test_outfile(ofname);
	file=fopen(ofname,"w");
	if (verbosity&VER_INPUT)
	  fprintf(stderr,"Opened %s for writing\n",ofname);
	if (stdo && (iter == iterations)) {
	  if (verbosity&VER_INPUT)
	    fprintf(stderr,"Writing to stdout\n");
	}
      }
      for (n=0;n<length;n++) {
	if (stdo) {
	  if (verbosity&VER_USR2) {
	    for (i=0;i<comp;i++) 
	      fprintf(stdout,"%e ",series[i][n]*interval[i]+min[i]);
	    fprintf(stdout,"%u\n",nmf[n]);
	  }
	  else {
	    fprintf(stdout,"%e",series[0][n]*interval[0]+min[0]);
	    for (i=1;i<comp;i++)
	      fprintf(stdout," %e",series[i][n]*interval[i]+min[i]);
	    fprintf(stdout,"\n");
	  }
	}
	if (!stdo || (verbosity&VER_USR1)) {
	  if (verbosity&VER_USR2) {
	    for (i=0;i<comp;i++) 
	      fprintf(file,"%e ",series[i][n]*interval[i]+min[i]);
	    fprintf(file,"%u\n",nmf[n]);
	  }
	  else {
	    fprintf(file,"%e",series[0][n]*interval[0]+min[0]);
	    for (i=1;i<comp;i++)
	      fprintf(file," %e",series[i][n]*interval[i]+min[i]);
	    fprintf(file,"\n");
	  }
	}
      }
      if (!stdo || (verbosity&VER_USR1))
	fclose(file);
    }
  }

  /*cleaning up */
  for (i=0;i<comp;i++) {
    free(series[i]);
    free(nf[i]);
    free(corr[i]);
  }
  free(series);
  free(nf);
  free(corr);

  for (i=0;i<2;i++)
    free(indexes[i]);
  free(indexes);

  free(list);
  free(nmf);
  free(interval);
  free(min);

  for (i=0;i<BOX;i++)
    free(box[i]);
  free(box);

  if (outfile != NULL)
    free(outfile);
  if (ofname != NULL)
    free(ofname);
  if (infile != NULL)
    free(infile);
  if (column != NULL)
    free(column);
  /* end cleaning up */

  return 0;
}
