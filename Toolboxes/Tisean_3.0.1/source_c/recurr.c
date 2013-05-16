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
/*Author: Rainer Hegger Last modified:  Sep 16, 2004 */
/* Sep 16, 2004: Change of index in output. before 0->N-1 now 1->N
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "This programs makes a recurrence plot for the data."

#define BOX 1024

unsigned long length=ULONG_MAX,exclude=0;
unsigned int embed=2,dim=1,delay=1;
unsigned int verbosity=0xff;
double eps=1.e-3,fraction=1.0;
char dimset=0;
char *columns;
char *outfile=NULL,stdo=1;
char *infile=NULL;
char epsset=0;

double **series;
long **box,*list;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"Usage: %s [options]\n",progname);
  fprintf(stderr,"Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to use [Default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [Default: 0]\n");
  fprintf(stderr,"\t-c columns to read [Default: 1]\n");
  fprintf(stderr,"\t-m # of components,embedding dimension [Default: 1,2]\n");
  fprintf(stderr,"\t-d delay [Default: 1]\n");
  fprintf(stderr,"\t-r size of the neighbourhood "
	  "[Default: (data interval)/1000]\n");
  fprintf(stderr,"\t-%% print only a percentage of points found [Default: "
	  " 100.0]\n");
  fprintf(stderr,"\t-o output file name [Default: 'datafile'.rec\n"
	  "\t\twithout -o: stdout]\n");
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
  if ((out=check_option(in,n,'c','s')) != NULL)
    columns=out;
  if ((out=check_option(in,n,'m','2')) != NULL) {
    sscanf(out,"%u,%u",&dim,&embed);
    dimset=1;
  }
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    epsset=1;
    sscanf(out,"%lf",&eps);
  }
  if ((out=check_option(in,n,'%','f')) != NULL) {
    sscanf(out,"%lf",&fraction);
    fraction /= 100.0;
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void lfind_neighbors(void)
{
  int i,i1,i2,j,j1,ke,ked,kd;
  int ibox=BOX-1;
  long n,element;
  double dx,epsinv;
  char toolarge;
  FILE *fout=NULL;

  epsinv=1./eps;
  rnd_init(0x9834725L);

  if (!stdo) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  for (n=(embed-1)*delay;n<length;n++) {
    i=(int)(series[0][n]*epsinv)&ibox;
    j=(int)(series[dim-1][n]*epsinv)&ibox;
    for (i1=i-1;i1<=i+1;i1++) {
      i2=i1&ibox;
      for (j1=j-1;j1<=j+1;j1++) {
	element=box[i2][j1&ibox];
	while (element > n) {
	  toolarge=0;
	  for (ke=0;ke<embed;ke++) {
	    ked=ke*delay;
	    for (kd=0;kd<dim;kd++) {
	      dx=fabs(series[kd][n-ked]-series[kd][element-ked]);
	      if (dx >= eps) {
		toolarge=1;
		break;
	      }
	    }
	    if (toolarge)
	      break;
	  }
	  if (!toolarge)
	    if (((double)rnd69069()/ULONG_MAX) <= fraction) {
	      if (!stdo)
		fprintf(fout,"%ld %ld\n",n+1,element+1);
	      else
		fprintf(stdout,"%ld %ld\n",n+1,element+1);
	    }
	  element=list[element];
	}
      }
    }
  }
  if (!stdo)
    fclose(fout);
}

int main(int argc,char **argv)
{
  long i;
  char stdi=0;
  double min,max,maxmax;

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
      strcat(outfile,".rec");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.rec");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (columns == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",dimset,
                                      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,columns,
                                      dimset,verbosity);

  maxmax=0.0;
  for (i=0;i<dim;i++) {
    rescale_data(series[i],length,&min,&max);
    if (max > maxmax)
      maxmax=max;
  }

  if (epsset)
    eps /= maxmax;

  check_alloc(list=(long*)malloc(sizeof(long)*length));
  check_alloc(box=(long**)malloc(sizeof(long*)*BOX));
  for (i=0;i<BOX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*BOX));

  make_multi_box(series,box,list,length,BOX,dim,embed,delay,eps);
  lfind_neighbors();

  return 0;
}
