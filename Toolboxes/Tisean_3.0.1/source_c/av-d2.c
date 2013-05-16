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
/*Author: Rainer Hegger Last modified: Sep 3, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Smoothes the output of the d2 program"

#define MAXLENGTH 1000

unsigned int maxdim=UINT_MAX,mindim=1;
unsigned int verbosity=0xff;
int aver=1;
char rescaled=0;
char stout=1;
char *outfile=NULL;
char *infile=NULL;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"Usage: %s [options]\n",progname);
  fprintf(stderr,"Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible datafile.\nStdin does NOT work.\n");
  fprintf(stderr,"\t-m dimension to start with [Default: 1]\n");
  fprintf(stderr,"\t-M dimension to end with [Default: whole file]\n");
  fprintf(stderr,"\t-a n average over (2n+1) values [Default: 1]\n");
  fprintf(stderr,"\t-E use rescaled data for the length scales\n\t\t"
	  "[Default: use units of data]\n");
  fprintf(stderr,"\t-o name of output file [Default: stdout,\n\t\t"
	  "-o without value means 'datafile'.av]\n");
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
    sscanf(out,"%u",&mindim);
  if ((out=check_option(in,n,'M','u')) != NULL)
    sscanf(out,"%u",&maxdim);
  if ((out=check_option(in,n,'a','u')) != NULL)
    sscanf(out,"%u",&aver);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'E','n')) != NULL)
    rescaled=1;
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char instr[1024];
  char *form1="%lf%lf",*form2="%*lf%lf%lf";
  char empty=0;
  unsigned int howmany,size=1;
  int j,k;
  long dim;
  double *eps,*y;
  double avy,aveps,norm;
  FILE *file,*fout=NULL;

  if ((argc < 2) || scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,0L,verbosity);
  if (infile == NULL) {
    fprintf(stderr,"You have to give a datafile. Exiting!\n");
    exit(127);
  }
  if (outfile == NULL) {
    check_alloc(outfile=(char*)calloc(strlen(infile)+4,(size_t)1));
    sprintf(outfile,"%s.av",infile);
  }
  
  check_alloc(eps=(double*)malloc(sizeof(double)*MAXLENGTH));
  check_alloc(y=(double*)malloc(sizeof(double)*MAXLENGTH));
  
  file=fopen(infile,"r");
  
  if (!stout) {
    test_outfile(outfile);
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  if (mindim > maxdim)
    mindim=maxdim;
  norm=2.0*aver+1.0;

  while (fgets(instr,1024,file) != NULL) {
    if (strlen(instr) != 1) {
      if (instr[0] == '#') {
	if (strstr(instr,"m= ") != NULL) {
	  sscanf(instr,"%*s %ld",&dim);
	  if ((dim >= mindim) && (dim <= maxdim)) {
	    howmany=0;
	    empty=0;
	    do {
	      if (fgets(instr,1024,file) == NULL)
		exit(127);
	      if (strlen(instr) == 1)
		empty=1;
	      if (!empty && (instr[0] != '#')) {
		if (!rescaled)
		  sscanf(instr,form1,&eps[howmany],&y[howmany]);
		else
		  sscanf(instr,form2,&y[howmany],&eps[howmany]);
		howmany++;
		if (!(howmany%MAXLENGTH)) {
		  check_alloc(realloc(eps,size*MAXLENGTH*sizeof(double)));
		  check_alloc(realloc(y,size*MAXLENGTH*sizeof(double)));
		  size++;
		}
	      }
	    } while (!empty);
	    for (k=aver;k<howmany-aver;k++) {
	      avy=aveps=0.0;
	      for (j= -aver;j<=aver;j++) {
		avy += y[k+j];
		aveps += eps[k+j];
	      }
	      if (!stout)
		fprintf(fout,"%e %e\n",aveps/norm,avy/norm);
	      else
		fprintf(stdout,"%e %e\n",aveps/norm,avy/norm);
	    }
	    if (!stout)
	      fprintf(fout,"\n");
	    else
	      fprintf(stdout,"\n");
	  }
	}
      }
    }
  }

  if (outfile != NULL)
    free(outfile);
  if (infile != NULL)
    free(infile);
  free(eps);
  free(y);

  return 0;
}
