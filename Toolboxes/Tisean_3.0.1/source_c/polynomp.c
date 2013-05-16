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
/*Author: Rainer Hegger. Last modified Sep 5, 2004 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "routines/tsa.h"

#define WID_STR "Fits a polynomial to the data."

char *outfile=NULL,stdo=1;
char *parin=NULL,*infile=NULL;
unsigned long length=ULONG_MAX,insample=ULONG_MAX,exclude=0;
unsigned long plength=UINT_MAX;
unsigned long step=1000;
unsigned int column=1,dim=2,delay=1,down_to=1;
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
  fprintf(stderr,"\t-L length of forecasted series [default: %lu]\n",step);
  fprintf(stderr,"\t-p name of parameter file [default: parameter.pol]\n");
  fprintf(stderr,"\t-o output file name [default: 'datafile'.pbf]\n");
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
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'n','u')) != NULL)
    sscanf(out,"%lu",&insample);
  if ((out=check_option(in,n,'L','u')) != NULL)
    sscanf(out,"%lu",&step);
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
  unsigned long n,hn;
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
  
  for (n=(dim-1)*delay;n<insample-1;n++) {
    hn=n+1;
    for (i=0;i<plength;i++) {
      vec[i] += series[hn]*(h=polynom(n,i));
      for (j=i;j<plength;j++)
	mat[i][j] += polynom(n,j)*h;
    }
  }
  for (i=0;i<plength;i++) {
    vec[i] /= (insample-(dim-1)*delay-1);
    for (j=i;j<plength;j++)
      mat[j][i]=(mat[i][j]/=(insample-(dim-1)*delay)-1);
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

  for (n=i0+(dim-1)*delay;n<i1-1;n++) {
    h=0.0;
    for (i=0;i<plength;i++)
      h += param[i]*polynom(n,i);
    error += (series[n+1]-h)*(series[n+1]-h);
  }

  return sqrt(error/(i1-i0-(dim-1)*delay-1));
}

void make_cast(FILE *fcast)
{
  int i,j,hi;
  unsigned int k;
  double casted;
  
  for (i=0;i<=(dim-1)*delay;i++)
    series[i]=series[length-(dim-1)*delay+i-1];

  hi=(dim-1)*delay;
  for (i=1;i<=step;i++) {
    casted=0.0;
    for (k=0;k<plength;k++)
      casted += param[k]*polynom((unsigned long)((dim-1)*delay),k);
    if (!stdo) {
      fprintf(fcast,"%e\n",casted);
      fflush(fcast);
    }
    else {
      fprintf(stdout,"%e\n",casted);
      fflush(stdout);
    }
    for (j=0;j<(dim-1)*delay;j++)
      series[j]=series[j+1];
    series[hi]=casted;
  }
}

int main(int argc,char **argv)
{
  int i,j;
  char stdi=0,oose=1;
  double **dummy,withalli,withallo;
  double av,varianz;
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
      sprintf(outfile,"%s.pbf",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.pbf");
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
    exit(POLYNOMP__WRONG_PARAMETER_FILE);
  }
  fclose(file);

  dummy=(double**)get_multi_series(parin,&plength,0LU,
				   &dim,"",(char)"1",verbosity);
  
  check_alloc(order=(unsigned int**)malloc(sizeof(int*)*plength));
  for (i=0;i<plength;i++) {
    check_alloc(order[i]=(unsigned int*)malloc(sizeof(int)*dim));
    for (j=0;j<dim;j++)
      order[i][j]=(unsigned int)dummy[j][i];
  }

  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  variance(series,length,&av,&varianz);

  if (insample >= length) {
    insample=length;
    oose=0;
  }

  check_alloc(param=(double*)malloc(sizeof(double)*plength));

  make_fit();
  withalli=forecast_error(0LU,insample);
  withallo=0.0;
  if (oose)
    withallo=forecast_error(insample+1,length);

  if (stdo) {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"#FCE: %e %e\n",withalli/varianz,withallo/varianz);
    for (i=0;i<plength;i++) {
      fprintf(stdout,"# ");
      for (j=0;j<dim;j++)
	fprintf(stdout,"%u ",order[i][j]);
      fprintf(stdout,"%e\n",param[i]);
    }
    fflush(stdout);
  }
  else {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(file,"#FCE: %e %e\n",withalli/varianz,withallo/varianz);
    for (i=0;i<plength;i++) {
      fprintf(file,"# ");
      for (j=0;j<dim;j++)
	fprintf(file,"%u ",order[i][j]);
      fprintf(file,"%e\n",param[i]);
    }
    fflush(file);
  }
  
  make_cast(file);

  if (!stdo)
    fclose(file);

  return 0;
}
