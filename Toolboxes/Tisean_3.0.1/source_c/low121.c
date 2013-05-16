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
/*Author: Rainer Hegger Last modified: Dec 17, 2001 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Simple lowpass filter in the time domain"


unsigned long length=ULONG_MAX,exclude=0;
unsigned int column=1,iterations=1;
unsigned int verbosity=0x1;
char *outfile=NULL,stdo=1;
char *infile=NULL;

double *series,*new;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"Usage: %s [options]\n",progname);
  fprintf(stderr,"Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of points to use [Default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [Default: 0]\n");
  fprintf(stderr,"\t-c column to read [Default: 1]\n");
  fprintf(stderr,"\t-i # of iterations [Default: 1]\n");
  fprintf(stderr,"\t-V verbosity level [Default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
          "2='+ print each iteration to a separate file\n");
  fprintf(stderr,"\t-o output file name(s) [Default: 'datafile'.low.n,\n\t\t"
	  "where n is the number of the iteration.\n\t\t"
	  "without -o the last iteration is written to stdout.]\n");
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
  if ((out=check_option(in,n,'i','u')) != NULL)
    sscanf(out,"%u",&iterations);
  if ((out=check_option(in,n,'V','d')) != NULL)
    sscanf(out,"%d",&verbosity);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  char *ofname;
  unsigned long i;
  unsigned int iter;
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
      check_alloc(ofname=(char*)calloc(strlen(infile)+9,(size_t)1));
      sprintf(outfile,"%s.low",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      check_alloc(ofname=(char*)calloc((size_t)14,(size_t)1));
      sprintf(outfile,"stdin.low");
    }
  }
  else
    check_alloc(ofname=(char*)calloc(strlen(outfile)+10,(size_t)1));
  
  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  check_alloc(new=(double*)malloc(sizeof(double)*length));
  
  if (verbosity&VER_USR1) {
    for (iter=1;iter<=iterations;iter++) {
      new[0]=(2.0*series[0]+2.0*series[1])/4.0;
      new[length-1]=(2.0*series[length-1]+2.0*series[length-2])/4.0;
      for (i=1;i<length-1;i++)
	new[i]=(series[i-1]+2.0*series[i]+series[i+1])/4.0;
	sprintf(ofname,"%s.%d",outfile,iter);
	test_outfile(ofname);
	file=fopen(ofname,"w");
	if (verbosity&VER_INPUT)
	  fprintf(stderr,"Opened %s for writing\n",ofname);
	if (stdo && (iter == iterations)) {
	  if (verbosity&VER_INPUT)
	    fprintf(stderr,"Writing to stdout\n");
	}
	for (i=0;i<length;i++) {
	  if (stdo && (iter == iterations))
	    fprintf(stdout,"%e\n",series[i]=new[i]);
	  fprintf(file,"%e\n",series[i]=new[i]);
	}
	fclose(file);
    }
  }
  else {
    for (iter=1;iter<=iterations;iter++) {
      new[0]=(2.0*series[0]+2.0*series[1])/4.0;
      new[length-1]=(2.0*series[length-1]+2.0*series[length-2])/4.0;
      for (i=1;i<length-1;i++)
	new[i]=(series[i-1]+2.0*series[i]+series[i+1])/4.0;
      for (i=0;i<length;i++)
	series[i]=new[i];
    }
    if (!stdo) {
      sprintf(ofname,"%s.%d",outfile,iterations);
      file=fopen(ofname,"w");
      if (verbosity&VER_INPUT)
	fprintf(stderr,"Opened %s for writing\n",ofname);
      for (i=0;i<length;i++)
	fprintf(file,"%e\n",series[i]);
      fclose(file);
    }
    else {
      if (verbosity&VER_INPUT)
	fprintf(stderr,"Writing to stdout\n");
      for (i=0;i<length;i++)
	fprintf(stdout,"%e\n",series[i]);
    }
  }

  return 0;
}
