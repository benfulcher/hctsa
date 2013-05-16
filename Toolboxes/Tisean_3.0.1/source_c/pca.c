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
/*Author: Rainer Hegger Last modified: Jul 26, 2004 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Performs a global PCA"

unsigned long LENGTH=ULONG_MAX,exclude=0;
unsigned int DIM=2,EMB=1,dimemb,LDIM=2,DELAY=1;
unsigned int verbosity=0xff;
char *outfile=NULL,stout=1,dim_set=0;
unsigned int what_to_write=0,write_values=1,write_vectors=0;
unsigned int write_comp=0,write_proj=0;
unsigned int projection_set=0;
char *infile=NULL,dimset=0,*column=NULL;
double **series;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to use [Default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignore [Default: 0]\n");
  fprintf(stderr,"\t-c columns to read [Default: 2]\n");
  fprintf(stderr,"\t-m columns,embedding dim. to use [Default: 2,1]\n");
  fprintf(stderr,"\t-d delay to use [Default: 1]\n");
  fprintf(stderr,"\t-q projection dimension [Default: no projection]\n");
  fprintf(stderr,"\t-W # what to write: [Default: 0]\n"
	  "\t\t0 write eigenvalues only\n"
	  "\t\t1 write eigenvectors\n"
	  "\t\t2 write (projected) pca components\n"
	  "\t\t3 write projected data\n");
  fprintf(stderr,"\t-o output file name \n\t\t[Default: stdout; -o without "
	  "value means 'datafile'.pca]\n");
  fprintf(stderr,"\t-V verbosity level [Default: 1]\n\t\t"
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
  if ((out=check_option(in,n,'c','s')) != NULL)
    column=out;
  if ((out=check_option(in,n,'m','2')) != NULL) {
    sscanf(out,"%u,%u",&DIM,&EMB);
    dimset=1;
  }
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(in,n,'q','u')) != NULL) {
    sscanf(out,"%u",&LDIM);
    projection_set=1;
  }
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'W','u')) != NULL) {
    sscanf(out,"%u",&what_to_write);
    switch(what_to_write) {
    case 0: write_values=1;break;
    case 1: write_values=0;write_vectors=1;break;
    case 2: write_values=0;write_comp=1;break;
    case 3: write_values=0;write_proj=1;break;
    default: {
      fprintf(stderr,"Wrong value for the -W flag. Exiting!\n");
      exit(127);
    }
    }
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void ordne(double *lyap,int *ord)
{
  long i,j,maxi;
  double max;
  
  for (i=0;i<dimemb;i++)
    ord[i]=i;

  for (i=0;i<dimemb-1;i++)
    for (j=i+1;j<dimemb;j++)
      if (lyap[i] < lyap[j]) {
	max=lyap[i];
	lyap[i]=lyap[j];
	lyap[j]=max;
	maxi=ord[i];
	ord[i]=ord[j];
	ord[j]=maxi;
      }
}

void make_pca(double *av)
{
  unsigned int i,j,k,i1,i2,j1,j2,k1,k2;
  int *ord;
  double **mat,*matarray,*eig,*sp,hsp=0.0;
  FILE *fout=NULL;

  check_alloc(ord=(int*)malloc(sizeof(int)*dimemb));
  check_alloc(eig=(double*)malloc(sizeof(double)*dimemb));
  check_alloc(matarray=(double*)malloc(sizeof(double)*dimemb*dimemb));
  check_alloc(mat=(double**)malloc(sizeof(double*)*dimemb));
  for (i=0;i<dimemb;i++)
    mat[i]=(double*)(matarray+i*dimemb);

  
  for (i=0;i<dimemb;i++) {
    i1=i/EMB;
    i2=(i%EMB)*DELAY;
    for (j=i;j<dimemb;j++) {
      j1=j/EMB;
      j2=(j%EMB)*DELAY;
      mat[i][j]=0.0;
      for (k=(EMB-1)*DELAY;k<LENGTH;k++)
	mat[i][j] += series[i1][k-i2]*series[j1][k-j2];
      mat[j][i]=(mat[i][j] /= (double)(LENGTH-(EMB-1)*DELAY));
    }
  }

  eigen(mat,(unsigned long)dimemb,eig);
  ordne(eig,ord);
  
  if (!stout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  for (i=0;i<dimemb;i++)
    if (write_values) {
      if (stout)
	fprintf(stdout,"%d %e\n",i,eig[i]);
      else
	fprintf(fout,"%d %e\n",i,eig[i]);
    }
    else {
      if (verbosity) {
	if (stout)
	  fprintf(stdout,"#%d %e\n",i,eig[i]);
	else
	  fprintf(fout,"#%d %e\n",i,eig[i]);
      }
    }
  if (write_vectors) {
    for (i=0;i<dimemb;i++) {
      for (j=0;j<dimemb;j++) {
	j1=ord[j];
	if (stout)
	  fprintf(stdout,"%e ",mat[i][j1]);
	else
	  fprintf(fout,"%e ",mat[i][j1]);
      }
      if (stout)
	fprintf(stdout,"\n");
      else
	fprintf(fout,"\n");
    }
  }

  if (write_comp) {
    for (i=(EMB-1)*DELAY;i<LENGTH;i++) {
      for (j=0;j<LDIM;j++) {
	j1=ord[j];
	hsp=0.0;
	for (k=0;k<dimemb;k++) {
	  k1=k/EMB;
	  k2=(k%EMB)*DELAY;
	  hsp += mat[k][j1]*(series[k1][i-k2]+av[k1]);
	}
	if (stout)
	  fprintf(stdout,"%e ",hsp);
	else
	  fprintf(fout,"%e ",hsp);
      }
      if (stout)
	fprintf(stdout,"\n");
      else
	fprintf(fout,"\n");
    }
  }

  if (write_proj) {
    check_alloc(sp=(double*)malloc(sizeof(double)*LDIM));
    for (i=0;i<(EMB-1)*DELAY;i++) {
      for (j=0;j<DIM;j++)
	if (stout)
	  fprintf(stdout,"%e ",series[j][i]+av[j]);
	else
	  fprintf(fout,"%e ",series[j][i]+av[j]);
      if (stout)
	fprintf(stdout,"\n");
      else
	fprintf(fout,"\n");
    }
    for (i=(EMB-1)*DELAY;i<LENGTH;i++) {
      for (j=0;j<LDIM;j++) {
	j1=ord[j];
	sp[j]=0.0;
	for (k=0;k<dimemb;k++) {
	  k1=k/EMB;
	  k2=(k%EMB)*DELAY;
	  sp[j] += mat[k][j1]*series[k1][i-k2];
	}
      }
      for (j=0;j<DIM;j++) {
	hsp=0.0;
	for (k=0;k<LDIM;k++) {
	  k1=ord[k];
	  hsp += mat[j*EMB][k1]*sp[k];
	}
	if (stout)
	  fprintf(stdout,"%e ",hsp+av[j]);
	else
	  fprintf(fout,"%e ",hsp+av[j]);
      }
      if (stout)
	fprintf(stdout,"\n");
      else
	fprintf(fout,"\n");
    }
    free(sp);
  }

  if (!stout)
    fclose(fout);
}

int main(int argc,char **argv)
{
  char stdi=0;
  unsigned int i,j;
  double rms,*av;

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
      strcat(outfile,".pca");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.pca");
    }
  }
  if (!stout)
    test_outfile(outfile);

  if (column == NULL)
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&DIM,"",dimset,
                                      verbosity);
  else
    series=(double**)get_multi_series(infile,&LENGTH,exclude,&DIM,column,
                                      dimset,verbosity);
  dimemb=DIM*EMB;
  if (!projection_set)
    LDIM=dimemb;
  else {
    if (LDIM < 1) LDIM=1;
    if (LDIM > dimemb) LDIM=dimemb;
  }

  check_alloc(av=(double*)malloc(sizeof(double)*DIM));
  for (j=0;j<DIM;j++) {
    av[j]=rms=0.0;
    variance(series[j],LENGTH,&av[j],&rms);
    for (i=0;i<LENGTH;i++)
      series[j][i] -= av[j];
  }
  make_pca(av);

  return 0;
}
