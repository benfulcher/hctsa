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
/*Author: Rainer Hegger. Last modified: Sep 4, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Creates a parameter file containing all terms\n\t\
for a polynomial"


char *outfile=NULL;
unsigned int dim=2,order=3;
unsigned int verbosity=0xff;
FILE *file=NULL;

void make_parameter(unsigned int*,unsigned int, unsigned int);

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [Options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"\t-m embedding dimension [Default: %u]\n",dim);
  fprintf(stderr,"\t-p order of the polynomial [Default: %u]\n",order);
  fprintf(stderr,"\t-o parameter file [Default: parameter.pol]\n");
  fprintf(stderr,"\t-V verbosity level [Default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **in)
{
  char *out;

  if ((out=check_option(in,n,'m','u')) != NULL)
    sscanf(out,"%u",&dim);
  if ((out=check_option(in,n,'p','u')) != NULL)
    sscanf(out,"%u",&order);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    if (strlen(out) > 0)
      outfile=out;
  }
}

void make_parameter(unsigned int *par,unsigned int d,unsigned int sum)
{
  int i,j;
  
  for (i=0;i<=order;i++) {
    sum += i;
    if (sum <= order) {
      par[d]=i;
      if (d == 0) {
	for (j=0;j<dim;j++)
	  fprintf(file,"%u ",par[j]);
	fprintf(file,"\n");
      }
      else
	make_parameter(par,d-1,sum);
    }
    sum -= i;
  }
  par[d]=0;
}

int main(int argc,char **argv)
{
  unsigned int i,*params;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  
  if (outfile == NULL) {
    check_alloc(outfile=(char*)calloc((size_t)14,(size_t)1));
    sprintf(outfile,"parameter.pol");
  }
  test_outfile(outfile);

  check_alloc(params=(unsigned int*)malloc(sizeof(unsigned int)*dim));
  for (i=0;i<dim;i++)
    params[i]=0;

  file=fopen(outfile,"w");
  if (verbosity&VER_INPUT)
    fprintf(stderr,"Opened %s for writing\n",outfile);
  make_parameter(params,dim-1,0);
  fclose(file);

  return 0;
}
