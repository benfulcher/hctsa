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
/*Author: Rainer Hegger Last modified: Dec 17, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Determines the maxima (minima) of a possibly multivariate\
 time series"


unsigned long length=ULONG_MAX,exclude=0;
char *column=NULL;
unsigned int verbosity=0xff;
unsigned int dim=1;
unsigned int which=1;
double mintime=0.0;
char dimset=0;
char maxima=1;
char stdo=1;
char *outfile=NULL;
char *infile=NULL;

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
  fprintf(stderr,"\t-m dimension (# of components) [Default: 1]\n");
  fprintf(stderr,"\t-c columns to read [Default: 1,...,# of components]\n");
  fprintf(stderr,"\t-w which component to maxi(mini)mize [Default: 1]\n");
  fprintf(stderr,"\t-z determine minima instead of maxima [Default: maxima]\n");
  fprintf(stderr,"\t-t minimal required time between two extrema "
	  "[Default: 0.0]\n");
  fprintf(stderr,"\t-o output file name [Default: 'datafile'.ext,"
	  " without -o: stdout]\n");
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
    sscanf(out,"%lu",&length);
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'m','u')) != NULL) {
    sscanf(out,"%u",&dim);
    dimset=1;
  }
  if ((out=check_option(in,n,'c','s')) != NULL)
    column=out;
  if ((out=check_option(in,n,'w','u')) != NULL)
    sscanf(out,"%u",&which);
  if ((out=check_option(in,n,'z','n')) != NULL)
    maxima=0;
  if ((out=check_option(in,n,'t','f')) != NULL)
    sscanf(out,"%lf",&mintime);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  unsigned long i,j;
  double **series;
  double x[3],a,b,c,lasttime,nexttime,time;
  FILE *fout=NULL;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  which--;
  if (which > (dim-1)) {
    fprintf(stderr,"The component to maxi(mini)mize has to be smaller or equal"
	    "to the number\nof components! Exiting\n");
    exit(EXTREMA_STRANGE_COMPONENT);
  }
  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      sprintf(outfile,"%s.ext",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.ext");
    }
  }

  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,column,
				      dimset,verbosity);
  
  if (!stdo) {
    test_outfile(outfile);    
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }
  
  lasttime=0.0;
  x[0]=series[which][0];
  x[1]=series[which][1];
  for (i=2;i<length;i++) {
    x[2]=series[which][i];
    if (maxima) {
      if ((x[1] >= x[0]) && (x[1] > x[2])) {
	a=x[1];
	b=(x[2]-x[0])/2.0;
	c=(x[2]-2.0*x[1]+x[0])/2.0;
	time= -b/2.0/c;
	nexttime=(double)i-1.0+time;
	if ((nexttime-lasttime) >= mintime) {
	  for (j=0;j<dim;j++) {
	    a=series[j][i-1];
	    b=(series[j][i]-series[j][i-2])/2.0;
	    c=(series[j][i]-2.0*series[j][i-1]+series[j][i-2])/2.0;
	    if (!stdo)
	      fprintf(fout,"%e ",a+b*time+c*sqr(time));
	    else
	      fprintf(stdout,"%e ",a+b*time+c*sqr(time));
	  }
	  if (!stdo)
	    fprintf(fout,"%e\n",nexttime-lasttime);
	  else
	    fprintf(stdout,"%e\n",nexttime-lasttime);
	  lasttime=nexttime;
	}
      }
    }
    else {
      if ((x[1] <= x[0]) && (x[1] < x[2])) {
	a=x[1];
	b=(x[2]-x[0])/2.0;
	c=(x[2]-2.0*x[1]+x[0])/2.0;
	time= -b/2.0/c;
	nexttime=(double)i-1.0+time;
	if ((nexttime-lasttime) >= mintime) {
	  for (j=0;j<dim;j++) {
	    a=series[j][i-1];
	    b=(series[j][i]-series[j][i-2])/2.0;
	    c=(series[j][i]-2.0*series[j][i-1]+series[j][i-2])/2.0;
	    if (!stdo)
	      fprintf(fout,"%e ",a+b*time+c*sqr(time));
	    else
	      fprintf(stdout,"%e ",a+b*time+c*sqr(time));
	  }
	  if (!stdo)
	    fprintf(fout,"%e\n",nexttime-lasttime);
	  else
	    fprintf(stdout,"%e\n",nexttime-lasttime);
	  lasttime=nexttime;
	}
      }
    }
    x[0]=x[1];
    x[1]=x[2];
  }
  if (!stdo)
    fclose(fout);

  if (infile != NULL)
    free(infile);
  if (outfile != NULL)
    free(outfile);
  for (i=0;i<dim;i++)
    free(series[i]);
  free(series);
  
  return 0;
}
