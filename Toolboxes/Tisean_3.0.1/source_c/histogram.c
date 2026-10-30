/*
 *   This file is part of TISEAN
 *
 *   Copyright (c) 1998-2015 Rainer Hegger, Holger Kantz, Thomas Schreiber
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
/*Author: Rainer Hegger. Last modified May 16, 2014*/

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Creates a histogram of a onedimensional dataset"

unsigned long length=ULONG_MAX;
unsigned long base=50;
unsigned long exclude=0;
unsigned int column=1;
unsigned int verbosity=0xff;
double size;
char my_stdout=1,gotsize=0,density=0;
char *outfile=NULL;
char *infile=NULL;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted as a"
	  " possible datafile.\nIf no datafile is given stdin is read. "
	  " Just - also means stdin\n");
  fprintf(stderr,"\t-l length of file [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default %ld]\n",exclude);
  fprintf(stderr,"\t-c column to read [default %d]\n",column);
  fprintf(stderr,"\t-b # of intervals [default %ld]\n",base);
  fprintf(stderr,"\t-D output densities not relative frequencies"
	  " [default not set]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.dat ;"
	  " If no -o is given: stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **str)
{
  char *out;

  if ((out=check_option(str,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(str,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(str,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(str,n,'b','u')) != NULL)
    sscanf(out,"%lu",&base);
  if ((out=check_option(str,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(str,n,'D','n')) != NULL)
    density=1;
  if ((out=check_option(str,n,'o','o')) != NULL) {
    my_stdout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  unsigned long i,j;
  double x,norm,size=1.0,size2=1.0;
  double min,max;
  double *series;
  double average,var;
  long *box;
  FILE *fout;

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
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,1));
      strcpy(outfile,infile);
      strcat(outfile,".his");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,1));
      strcpy(outfile,"stdin.his");
    }
  }
  if (!my_stdout)
    test_outfile(outfile);

  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  variance(series,length,&average,&var);

  min=max=series[0];
  for (i=1;i<length;i++) {
    if (series[i] < min) min=series[i];
    else if (series[i] > max) max=series[i];
  }
  max -= min;

  for (i=0;i<length;i++)
    series[i]=(series[i]-min);

  if (base > 0) {
    check_alloc(box=(long*)malloc(sizeof(long)*base));
    for (i=0;i<base;i++)
      box[i]=0;
    size=1./base;
    size2=(1.0-size/2.0)*max;
    for (i=0;i<length;i++) {
      if (series[i] > size2)
	series[i]=size2;
      j=(long)(series[i]*base/max);
      box[j]++;
    }
  }

  if (!density)
    norm=1.0/(double)length;
  else
    norm=1.0/(double)length*(double)base/max;

  if (!my_stdout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(fout,"#interval of data: [%e:%e]\n",min,max+min);
    fprintf(fout,"#average= %e\n",average);
    fprintf(fout,"#standard deviation= %e\n",var);
    for (i=0;i<base;i++) {
      x=(double)(i*size);
      fprintf(fout,"%e %e\n",(x+size/2.0)*max+min,(double)box[i]*norm);
    }
    fclose(fout);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"#interval of data: [%e:%e]\n",min,max+min);
    fprintf(stdout,"#average= %e\n",average);
    fprintf(stdout,"#standard deviation= %e\n",var);
    for (i=0;i<base;i++) {
      x=(double)(i*size);
      fprintf(stdout,"%e %e\n",(x+size/2.0)*max+min,(double)box[i]*norm);
      fflush(stdout);
    }
  }
  return 0;
}
