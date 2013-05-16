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
/*Author: Rainer Hegger. Last modified: Dec 10, 2005 */
/*Changes:
  12/10/05: It's multivariate now
  12/16/05: Scaled <eps> and sigma(eps)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "routines/tsa.h"

#define WID_STR "Determines the fraction of false nearest neighbors."

char *outfile=NULL;
char *infile=NULL;
char stdo=1,dimset=0;
char *column=NULL;
unsigned long length=ULONG_MAX,exclude=0,theiler=0;
unsigned int delay=1,maxdim=5,minemb=1;
unsigned int comp=1,maxemb=5;
unsigned int verbosity=0xff;
double rt=2.0;
double eps0=1.0e-5;
double **series;
double aveps,vareps;
double varianz;

#define BOX 1024
int ibox=BOX-1;
long **box,*list;
unsigned int *vcomp,*vemb;
unsigned long toolarge;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data [default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default: 0]\n");
  fprintf(stderr,"\t-c columns to read [default: 1]\n");
  fprintf(stderr,"\t-m min. test embedding dimension [default: %u]\n",minemb);
  fprintf(stderr,"\t-M # of components,max. emb. dim. [default: %u,%u]\n",
	  comp,maxemb);
  fprintf(stderr,"\t-d delay [default: 1]\n");
  fprintf(stderr,"\t-f escape factor [default: %.2lf]\n",rt);
  fprintf(stderr,"\t-t theiler window [default: 0]\n");
  fprintf(stderr,"\t-o output file [default: 'datafile'.fnn; without -o"
	  " stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default: 3]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
          "2='+ information about the current state\n");
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
    column=out;
  if ((out=check_option(in,n,'m','u')) != NULL)
    sscanf(out,"%u",&minemb);
  if ((out=check_option(in,n,'M','2')) != NULL) {
    sscanf(out,"%u,%u",&comp,&maxemb);
    maxdim=comp*(maxemb+1);
    dimset=1;
  }
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(in,n,'f','f')) != NULL)
    sscanf(out,"%lf",&rt);
  if ((out=check_option(in,n,'t','u')) != NULL)
    sscanf(out,"%lu",&theiler);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void mmb(unsigned int hdim,unsigned int hemb,double eps)
{
  unsigned long i;
  long x,y;

  for (x=0;x<BOX;x++)
    for (y=0;y<BOX;y++)
      box[x][y] = -1;

  for (i=0;i<length-(maxemb+1)*delay;i++) {
    x=(long)(series[0][i]/eps)&ibox;
    y=(long)(series[hdim][i+hemb]/eps)&ibox;
    list[i]=box[x][y];
    box[x][y]=i;
  }
}

char find_nearest(long n,unsigned int dim,double eps)
{
  long x,y,x1,x2,y1,i,i1,ic,ie;
  long element,which= -1;
  double dx,maxdx,mindx=1.1,hfactor,factor;

  ic=vcomp[dim];
  ie=vemb[dim];
  x=(long)(series[0][n]/eps)&ibox;
  y=(long)(series[ic][n+ie]/eps)&ibox;
  
  for (x1=x-1;x1<=x+1;x1++) {
    x2=x1&ibox;
    for (y1=y-1;y1<=y+1;y1++) {
      element=box[x2][y1&ibox];
      while (element != -1) {
	if (labs(element-n) > theiler) {
	  maxdx=fabs(series[0][n]-series[0][element]);
	  for (i=1;i<=dim;i++) {
	    ic=vcomp[i];
	    i1=vemb[i];
	    dx=fabs(series[ic][n+i1]-series[ic][element+i1]);
	    if (dx > maxdx)
	      maxdx=dx;
	  }
	  if ((maxdx < mindx) && (maxdx > 0.0)) {
	    which=element;
	    mindx=maxdx;
	  }
	}
	element=list[element];
      }
    }
  }

  if ((which != -1) && (mindx <= eps) && (mindx <= varianz/rt)) {
    aveps += mindx;
    vareps += mindx*mindx;
    factor=0.0;
    for (i=1;i<=comp;i++) {
      ic=vcomp[dim+i];
      ie=vemb[dim+i];
      hfactor=fabs(series[ic][n+ie]-series[ic][which+ie])/mindx;
      if (hfactor > factor) 
	factor=hfactor;
    }
    if (factor > rt)
      toolarge++;
    return 1;
  }
  return 0;
}

int main(int argc,char **argv)
{
  char stdi=0;
  FILE *file=NULL;
  double min,inter=0.0,ind_inter,epsilon,av,ind_var;
  char *nearest,alldone;
  long i;
  unsigned int dim,emb;
  unsigned long donesofar;

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
      strcat(outfile,".fnn");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.fnn");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&comp,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&comp,column,
				      dimset,verbosity);

  for (i=0;i<comp;i++) {
    rescale_data(series[i],length,&min,&ind_inter);
    variance(series[i],length,&av,&ind_var);
    if (i == 0) {
      varianz=ind_var;
      inter=ind_inter;
    }
    else {
      varianz=(varianz>ind_var)?ind_var:varianz;
      inter=(inter<ind_inter)?ind_inter:inter;
    }
  }

  check_alloc(list=(long*)malloc(sizeof(long)*length));
  check_alloc(nearest=(char*)malloc(length));
  check_alloc(box=(long**)malloc(sizeof(long*)*BOX));
  for (i=0;i<BOX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*BOX));

  if (!stdo) {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }
  check_alloc(vcomp=(unsigned int*)malloc(sizeof(int)*(maxdim)));
  check_alloc(vemb=(unsigned int*)malloc(sizeof(int)*(maxdim)));
  for (i=0;i<maxdim;i++) {
    if (comp == 1) {
      vcomp[i]=0;
      vemb[i]=i;
    }
    else {
      vcomp[i]=i%comp;
      vemb[i]=(i/comp)*delay;
    }
  }
  for (emb=minemb;emb<=maxemb;emb++) {
    dim=emb*comp-1;
    epsilon=eps0;
    toolarge=0;
    alldone=0;
    donesofar=0;
    aveps=0.0;
    vareps=0.0;
    for (i=0;i<length;i++)
      nearest[i]=0;
    if (verbosity&VER_USR1)
      fprintf(stderr,"Start for dimension=%u\n",dim+1);
    while (!alldone && (epsilon < 2.*varianz/rt)) {
      alldone=1;
      mmb(vcomp[dim],vemb[dim],epsilon);
      for (i=0;i<length-maxemb*delay;i++)
	if (!nearest[i]) {
	  nearest[i]=find_nearest(i,dim,epsilon);
	  alldone &= nearest[i];
	  donesofar += (unsigned long)nearest[i];
	}
      if (verbosity&VER_USR1)
	fprintf(stderr,"Found %lu up to epsilon=%e\n",donesofar,epsilon*inter);
      epsilon*=sqrt(2.0);
      if (!donesofar)
	eps0=epsilon;
    }
    if (donesofar == 0) {
      fprintf(stderr,"Not enough points found!\n");
      exit(FALSE_NEAREST_NOT_ENOUGH_POINTS);
    }
    aveps *= (1./(double)donesofar);
    vareps *= (1./(double)donesofar);
    if (stdo) {
      fprintf(stdout,"%u %e %e %e\n",dim+1,(double)toolarge/(double)donesofar,
	      aveps*inter,sqrt(vareps)*inter);
      fflush(stdout);
    }
    else {
      fprintf(file,"%u %e %e %e\n",dim+1,(double)toolarge/(double)donesofar,
	      aveps*inter,sqrt(vareps)*inter);
      fflush(file);
    }
  }
  if (!stdo)
    fclose(file);

  if (infile != NULL)
    free(infile);
  if (outfile != NULL)
    free(outfile);
  free(series);
  free(list);
  free(nearest);
  for (i=0;i<BOX;i++)
    free(box[i]);
  free(box);

  return 0;
}
