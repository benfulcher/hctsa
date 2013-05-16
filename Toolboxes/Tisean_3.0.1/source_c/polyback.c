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
/*Author: Rainer Hegger. Last modified Sep 4, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "routines/tsa.h"

#define WID_STR "Does a backward elimination for a polynomial"

char *outfile=NULL,stdo=1;
char *parin=NULL,*infile=NULL;
unsigned long length=ULONG_MAX,insample=ULONG_MAX,exclude=0;
unsigned int plength=UINT_MAX;
unsigned int column=1,dim=2,delay=1,down_to=1,step=1;
unsigned int **order;
unsigned int verbosity=0xff;
double *series,*param;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"Usage: %s [Options]\n",progname);
  fprintf(stderr,"Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to use [default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default: %lu]\n",exclude);
  fprintf(stderr,"\t-c column to read [default: %u]\n",column);
  fprintf(stderr,"\t-m embedding dimension [default: %u]\n",dim);
  fprintf(stderr,"\t-d delay [default: %u]\n",delay);
  fprintf(stderr,"\t-n insample data [default: all]\n");
  fprintf(stderr,"\t-s steps to forecast [default: %u]\n",step);
  fprintf(stderr,"\t-# reduce down to # terms [default: %u]\n",down_to);
  fprintf(stderr,"\t-p name of parameter file [default: parameter.pol]\n");
  fprintf(stderr,"\t-o output file name [default: 'datafile'.pbe]\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
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
  if ((out=check_option(in,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(in,n,'m','u')) != NULL)
    sscanf(out,"%u",&dim);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(in,n,'n','u')) != NULL)
    sscanf(out,"%lu",&insample);
  if ((out=check_option(in,n,'#','u')) != NULL)
    sscanf(out,"%u",&down_to);
  if ((out=check_option(in,n,'s','u')) != NULL)
    sscanf(out,"%u",&step);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'p','s')) != NULL)
    parin=out;
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

double polynom(unsigned long act,unsigned int which)
{
  unsigned int i,j;
  double ret=1.0,h;
  
  for (i=0;i<dim;i++) {
    h=series[act-i*delay];
    for (j=0;j<order[which][i];j++)
      ret *= h;
  }
  
  return ret;
}

void make_fit(void)
{
  double **mat,*vec;
  double h;
  unsigned long n;
  unsigned int i,j;

  check_alloc(vec=(double*)malloc(sizeof(double)*plength));
  check_alloc(mat=(double**)malloc(sizeof(double*)*plength));
  for (i=0;i<plength;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*plength));

  for (i=0;i<plength;i++) {
    vec[i]=0.0;
    for (j=0;j<plength;j++)
      mat[i][j]=0.0;
  }
  
  for (n=(dim-1)*delay;n<insample-step;n++) {
    for (i=0;i<plength;i++) {
      vec[i] += series[n+step]*(h=polynom(n,i));
      for (j=i;j<plength;j++)
	mat[i][j] += polynom(n,j)*h;
    }
  }
  for (i=0;i<plength;i++) {
    vec[i] /= (insample-step-(dim-1)*delay);
    for (j=i;j<plength;j++)
      mat[j][i]=(mat[i][j]/=(insample-step-(dim-1)*delay));
  }
  
  solvele(mat,vec,plength);

  for (i=0;i<plength;i++)
    param[i]=vec[i];

  free(vec);
  for (i=0;i<plength;i++)
    free(mat[i]);
  free(mat);
}

double forecast_error(unsigned long i0,unsigned long i1)
{
  unsigned int i;
  unsigned long n;
  double h,error=0.0;

  for (n=i0+(dim-1)*delay;n<i1-step;n++) {
    h=0.0;
    for (i=0;i<plength;i++)
      h += param[i]*polynom(n,i);
    error += (series[n+step]-h)*(series[n+step]-h);
  }
  
  return sqrt(error/(i1-i0-step-(dim-1)*delay));
}

int main(int argc,char **argv)
{
  int i,j,k,l,hl,ibest,counter;
  char stdi=0,out_set=1,*parout;
  double **dummy,besti,besto,withalli,withallo,errori=0.,erroro=0.;
  double av,varianz;
  unsigned long hlength=ULONG_MAX;
  unsigned int **ini_params,*isout,offset;
  FILE *file,*fpars;

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
      sprintf(outfile,"%s.pbe",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.pbe");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (parin == NULL) {
    check_alloc(parin=(char*)calloc((size_t)14,(size_t)1));
    sprintf(parin,"parameter.pol");
  }
  file=fopen(parin,"r");
  if (file == NULL) {
    fprintf(stderr,"File %s does not exist. Exiting!\n",parin);
    exit(POLYBACK__WRONG_PARAMETER_FILE);
  }
  fclose(file);

  if (verbosity&VER_INPUT)
    fprintf(stderr,"Using %s as the parameter file\n",parin);
  dummy=(double**)get_multi_series(parin,&hlength,0LU,&dim,"",(char)1,
				   verbosity);

  offset=(unsigned int)(log((double)hlength)/log(10.0)+1.0);
  check_alloc(parout=(char*)calloc(strlen(parin)+offset+2,(size_t)1));
  
  check_alloc(ini_params=(unsigned int**)malloc(sizeof(int*)*hlength));
  for (i=0;i<hlength;i++) {
    check_alloc(ini_params[i]=(unsigned int*)malloc(sizeof(int)*dim));
    for (j=0;j<dim;j++)
      ini_params[i][j]=(unsigned int)dummy[j][i];
  }
  check_alloc(isout=(unsigned int*)malloc(sizeof(int)*hlength));

  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  variance(series,length,&av,&varianz);

  if (insample >= length) {
    insample=length;
    out_set=0;
  }

  check_alloc(order=(unsigned int**)malloc(sizeof(int*)*hlength));
  check_alloc(param=(double*)malloc(sizeof(double)*hlength));
  for (i=0;i<hlength;i++) {
    isout[i]=0;
    check_alloc(order[i]=(unsigned int*)malloc(sizeof(int)*dim));
    for (j=0;j<dim;j++)
      order[i][j]=ini_params[i][j];
  }
  plength=hlength;

  make_fit();
  withalli=forecast_error(0LU,insample);
  withallo=0.0;
  if (out_set)
    withallo=forecast_error(insample+1,length);

  if (stdo) {
    fprintf(stdout,"%lu %e %e\n",hlength,withalli/varianz,withallo/varianz);
    fflush(stdout);
  }
  else {
    file=fopen(outfile,"w");
    fprintf(file,"%lu %e %e\n",hlength,withalli/varianz,withallo/varianz);
    fflush(file);
  }
  free(param);
  for (i=0;i<plength;i++)
    free(order[i]);
  free(order);
  
  if ((down_to < 1) || (down_to > hlength))
    down_to=1;

  for (i=1;i<=hlength-down_to;i++) {
    plength=hlength-i;
    besti=besto=0.0;
    ibest= -1;
    check_alloc(order=(unsigned int**)malloc(sizeof(int*)*plength));
    check_alloc(param=(double*)malloc(sizeof(double)*plength));
    for (j=0;j<plength;j++) {
      check_alloc(order[j]=(unsigned int*)malloc(sizeof(int)*dim));
    }
    counter=plength;
    for (j=0;j<hlength;j++)
      if (!isout[j]) {
	isout[j]++;
	hl=0;
	for (k=0;k<hlength;k++) {
	  if (!isout[k]) {
	    for (l=0;l<dim;l++)
	      order[hl][l]=ini_params[k][l];
	    hl++;
	  }
	}
	make_fit();
	errori=forecast_error(0LU,insample);
	if (out_set)
	  erroro=forecast_error(insample+1,length);
	if (ibest == -1) {
	  besti=errori;
	  if (out_set)
	    besto=erroro;
	  ibest=j;
	}
	else {
	  if (out_set) {
	    if (erroro < besto) {
	      besto=erroro;
	      besti=errori;
	      ibest=j;
	    }
	  }
	  else {
	    if (errori < besti) {
	      besti=errori;
	      besto=erroro;
	      ibest=j;
	    }
	  }
	}
	isout[j]--;
      }
    isout[ibest]++;
    free(param);
    for (j=0;j<plength;j++)
      free(order[j]);
    free(order);
    if (stdo) {
      fprintf(stdout,"%u %e %e ",plength,besti/varianz,besto/varianz);
      for (j=0;j<dim;j++)
	fprintf(stdout,"%u ",ini_params[ibest][j]);
      fprintf(stdout,"\n");
      fflush(stdout);
    }
    else {
      fprintf(file,"%u %e %e ",plength,besti/varianz,besto/varianz);
      for (j=0;j<dim;j++)
	fprintf(file,"%u ",ini_params[ibest][j]);
      fprintf(file,"\n");
      fflush(file);
    }
    sprintf(parout,"%s.%u",parin,plength);
    fpars=fopen(parout,"w");
    for (j=0;j<hlength;j++)
      if (!isout[j]) {
	for (k=0;k<dim;k++)
	  fprintf(fpars,"%u ",ini_params[j][k]);
	fprintf(fpars,"\n");
      }
    fclose(fpars);
  }
 
  if (!stdo)
    fclose(file);
  return 0;
}
