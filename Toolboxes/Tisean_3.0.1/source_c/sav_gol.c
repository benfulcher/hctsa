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
/*Author: Rainer Hegger. Last modified May 27, 2000 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Savitzky-Golay filter: Filters the data or estimates\n\t\
filtered derivatives, respectively."

unsigned long length=ULONG_MAX,exclude=0;
unsigned int dim=1;
char dimset=0;
char *columns=NULL;
unsigned int nf=2,nb=2,power=2,deriv=0;
char *infile=NULL,*outfile=NULL,stdo=1;
unsigned int verbosity=(VER_INPUT|VER_FIRST_LINE);

double **series;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"  Usage: %s [options]\n",progname);
  fprintf(stderr,"  Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l datapoints [default is whole file]\n");
  fprintf(stderr,"\t-x exclude # points [default %ld]\n",exclude);
  fprintf(stderr,"\t-c columns [default 1]\n");
  fprintf(stderr,"\t-m no. of components [default %d]\n",dim);
  fprintf(stderr,"\t-n nb,nf [default %u,%u]\n",nb,nf);
  fprintf(stderr,"\t-p power of the polynomial [default %u]\n",power);
  fprintf(stderr,"\t-D order of the estimated derivative [default %u]\n",deriv);
  fprintf(stderr," \t-o outfile [default 'datafile'.sg; Without -o data"
	  " is written to stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
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
  if ((out=check_option(argv,n,'c','s')) != NULL)
    columns=out;
  if ((out=check_option(argv,n,'m','u')) != NULL) {
    sscanf(out,"%u",&dim);
    dimset=1;
  }
  if ((out=check_option(argv,n,'n','2')) != NULL)
    sscanf(out,"%u,%u",&nb,&nf);
  if ((out=check_option(argv,n,'p','u')) != NULL)
    sscanf(out,"%u",&power);
  if ((out=check_option(argv,n,'D','u')) != NULL)
    sscanf(out,"%u",&deriv);
  if ((out=check_option(argv,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

double** make_coeff(void)
{
  long i,j,k;
  double **mat,**imat,**rmat;
  
  check_alloc(mat=(double**)malloc(sizeof(double*)*(power+1)));
  for (i=0;i<=power;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*(power+1)));
  check_alloc(rmat=(double**)malloc(sizeof(double*)*(power+1)));
  for (i=0;i<=power;i++)
    check_alloc(rmat[i]=(double*)malloc(sizeof(double)*(nb+nf+1)));
  
  for (i=0;i<=power;i++)
    for (j=0;j<=power;j++) {
      mat[i][j]=0.0;
      for (k= -(int)nb;k<=(int)nf;k++)
	mat[i][j] += pow((double)k,(double)(i+j));
    }

  imat=invert_matrix(mat,(power+1));
  
  for (i=0;i<=power;i++)
    for (j=0;j<=(nb+nf);j++) {
      rmat[i][j]=0.0;
      for (k=0;k<=power;k++)
	rmat[i][j] += imat[i][k]*pow((double)(j-(int)nb),(double)k);
    }
  
  for (i=0;i<=power;i++) {
    free(mat[i]);
    free(imat[i]);
  }
  free(mat);
  free(imat);

  return rmat;
}

double make_norm(void)
{
  double ret=1.0;
  long i;

  for (i=2;i<=deriv;i++)
    ret *= (double)i;

  return 1.0/ret;
}

int main(int argc,char **argv)
{
  char stdi=0;
  long i,j,d;
  double **coeff,help,norm;
  FILE *fout;

  if (scan_help(argc,argv)) 
    show_options(argv[0]);
  
  scan_options(argc,argv);

  if (power >= (nb+nf+1)) {
    fprintf(stderr,"With these settings for the -n and -p flags,\nthe"
	    " system is underdetermined. Exiting\n\n");
    exit(SAV_GOL__UNDERDETERMINED);
  }
  if (deriv > power) {
    fprintf(stderr,"The order of the derivative must not be larger\nthan"
	    " the power of polynomial. Exiting\n\n");
    exit(SAV_GOL__TOO_LARGE_DERIVATIVE);
  }

#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+4,(size_t)1));
      sprintf(outfile,"%s.sg",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)9,(size_t)1));
      sprintf(outfile,"stdin.sg");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (columns == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",
				      dimset,verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,
				      columns,dimset,verbosity);
  
  coeff=make_coeff();
  norm=make_norm();

  if (stdo) {
    for (i=0;i<nb;i++) {
      for (d=0;d<dim;d++)
	fprintf(stdout,"%e ",(deriv==0)?series[d][i]:0.0);
      fprintf(stdout,"\n");
    }
    for (i=(long)nb;i<length-(long)nf;i++) {
      for (d=0;d<dim;d++) {
	help=0.0;
	for (j= -(long)nb;j<=(long)nf;j++)
	  help += coeff[deriv][j+nb]*series[d][i+j];
	fprintf(stdout,"%e ",help*norm);
      }
      fprintf(stdout,"\n");
    }
    for (i=length-(long)nf;i<length;i++) {
      for (d=0;d<dim;d++)
	fprintf(stdout,"%e ",(deriv==0)?series[d][i]:0.0);
      fprintf(stdout,"\n");
    }
  }
  else {
    fout=fopen(outfile,"w");
    for (i=0;i<nb;i++) {
      for (d=0;d<dim;d++)
	fprintf(fout,"%e ",(deriv==0)?series[d][i]:0.0);
      fprintf(fout,"\n");
    }
    for (i=(long)nb;i<length-(long)nf;i++) {
      for (d=0;d<dim;d++) {
	help=0.0;
	for (j= -(long)nb;j<=(long)nf;j++)
	  help += coeff[deriv][j+nb]*series[d][i+j];
	fprintf(fout,"%e ",help*norm);
      }
      fprintf(fout,"\n");
    }
    for (i=length-(long)nf;i<length;i++) {
      for (d=0;d<dim;d++)
	fprintf(fout,"%e ",(deriv==0)?series[d][i]:0.0);
      fprintf(fout,"\n");
    }
    fclose(fout);
  }

  for (i=0;i<dim;i++)
    free(series[i]);
  free(series);
  free(outfile);
  if (!stdi)
    free(infile);
  for (i=0;i<=power;i++)
    free(coeff[i]);
  free(coeff);

  return 0;
}
