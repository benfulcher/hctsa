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
/*Author: Rainer Hegger. Last modified: Dec 10, 2005 */
/*Changes:
  12/10/05: It's multivariate now
  12/16/05: Scaled <eps> and sigma(eps)
  03/08/09: delay was missing in delay embedding of univariate case
  hctsa: a fifth output column, the number of points that found a neighbor
         within the search range at that dimension (the denominator of the
         fraction of false nearest neighbors)
  hctsa: distances are compared with a tolerance (DIST_TOL, in the units of
         the data rescaled to [0,1]), so that the result does not depend on
         rounding: distances within DIST_TOL of each other are equal, and a
         pair closer than DIST_TOL is a repeated point (no neighbor). When
         several neighbors are equally near, each point contributes the
         fraction of them that are false (the expected outcome of choosing
         one of them at random), so the result does not depend on the order
         in which the neighbors are searched.
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
unsigned int delay=1,maxdim=6,minemb=1;
unsigned int comp=1,maxemb=5;
unsigned int verbosity=0xff;
double rt=2.0;
/* distances closer than this (data rescaled to [0,1]) are equal: some 10^4 times
   the rounding error of the rescaled data (about 10^-16), and far finer than any
   measured data resolve */
#define DIST_TOL 1.0e-12
double eps0=1.0e-5;
double **series;
double aveps,vareps;
double varianz;

#define BOX 1024
int ibox=BOX-1;
long **box,*list;
unsigned int *vcomp,*vemb;
double toolarge;

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

/* the max-norm distance between points n and element in dim+1 dimensions */
static double maxnorm(long n,long element,unsigned int dim)
{
  long i,ic,i1;
  double dx,maxdx;

  maxdx=fabs(series[0][n]-series[0][element]);
  for (i=1;i<=dim;i++) {
    ic=vcomp[i];
    i1=vemb[i];
    dx=fabs(series[ic][n+i1]-series[ic][element+i1]);
    if (dx > maxdx)
      maxdx=dx;
  }
  return maxdx;
}

char find_nearest(long n,unsigned int dim,double eps)
{
  long x,y,x1,x2,y1,i,ic,ie;
  long element,which= -1;
  unsigned long nnear,nfalse;
  double maxdx,mindx=1.1,dx,factor;

  ic=vcomp[dim];
  ie=vemb[dim];
  x=(long)(series[0][n]/eps)&ibox;
  y=(long)(series[ic][n+ie]/eps)&ibox;

  /* the distance to the nearest neighbor (repeated points excluded) */
  for (x1=x-1;x1<=x+1;x1++) {
    x2=x1&ibox;
    for (y1=y-1;y1<=y+1;y1++) {
      element=box[x2][y1&ibox];
      while (element != -1) {
	if (labs(element-n) > theiler) {
	  maxdx=maxnorm(n,element,dim);
	  if ((maxdx < mindx) && (maxdx > DIST_TOL)) {
	    which=element;
	    mindx=maxdx;
	  }
	}
	element=list[element];
      }
    }
  }

  /* every point within mindx+DIST_TOL is in the 3x3 boxes once mindx+DIST_TOL <= eps */
  if ((which != -1) && (mindx+DIST_TOL <= eps) && (mindx <= varianz/rt)) {
    aveps += mindx;
    vareps += mindx*mindx;
    /* the nearest neighbors: all within DIST_TOL of the nearest distance */
    nnear=nfalse=0;
    for (x1=x-1;x1<=x+1;x1++) {
      x2=x1&ibox;
      for (y1=y-1;y1<=y+1;y1++) {
	element=box[x2][y1&ibox];
	while (element != -1) {
	  if (labs(element-n) > theiler) {
	    maxdx=maxnorm(n,element,dim);
	    if ((maxdx > DIST_TOL) && (maxdx <= mindx+DIST_TOL)) {
	      factor=0.0;
	      for (i=1;i<=comp;i++) {
		ic=vcomp[dim+i];
		ie=vemb[dim+i];
		dx=fabs(series[ic][n+ie]-series[ic][element+ie]);
		if (dx > factor)
		  factor=dx;
	      }
	      nnear++;
	      /* false if the next coordinate separates them by more than rt times
		 their distance */
	      if (factor > rt*maxdx+DIST_TOL)
		nfalse++;
	    }
	  }
	  element=list[element];
	}
      }
    }
    toolarge += (double)nfalse/(double)nnear;
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

  /* Guard against unsigned underflow below (length-(maxemb+1)*delay and
     length-maxemb*delay, both computed as unsigned arithmetic): if the
     requested embedding dimension times the delay is not comfortably
     smaller than the data length, those subtractions wrap around to a
     huge positive value instead of going negative, and the main loops
     then index nearest[]/series[][] far out of bounds. Fail cleanly
     instead. */
  if ((unsigned long)(maxemb+1)*delay >= length) {
    fprintf(stderr,"false_nearest: maximal embedding dimension (%u) times "
	    "delay (%u) is too large for the data length (%lu). Reduce -M "
	    "or -d, or supply more data.\n",maxemb,delay,length);
    exit(FALSE_NEAREST_NOT_ENOUGH_POINTS);
  }

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
  /* Embedding coordinate i is component i%comp at lag (i/comp)*delay. For a
     scalar series (comp == 1) this is a delay embedding with lag i*delay.
     (TISEAN 3.0.0/3.0.1 as released special-cased comp == 1 with vemb[i]=i,
     so the -d option never set the lag between coordinates in the scalar
     case, only the loop bounds above; the authors fixed this in 2009 in an
     unreleased version, see the 03/08/09 note above.) */
  for (i=0;i<maxdim;i++) {
    vcomp[i]=i%comp;
    vemb[i]=(i/comp)*delay;
  }
  for (emb=minemb;emb<=maxemb;emb++) {
    dim=emb*comp-1;
    epsilon=eps0;
    toolarge=0.0;
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
      fprintf(stdout,"%u %e %e %e %lu\n",dim+1,toolarge/(double)donesofar,
	      aveps*inter,sqrt(vareps)*inter,donesofar);
      fflush(stdout);
    }
    else {
      fprintf(file,"%u %e %e %e %lu\n",dim+1,toolarge/(double)donesofar,
	      aveps*inter,sqrt(vareps)*inter,donesofar);
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
