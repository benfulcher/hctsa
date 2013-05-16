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
/*Author: Rainer Hegger Last modified: Jun 10, 2006 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Multivariate noise reduction using the GHKSS algorithm"


#define BOX (unsigned int)1024

unsigned long length=ULONG_MAX,exclude=0;
unsigned int dim,qdim=2,delay=1,minn=50,iterations=1,comp=1,embed=5;
unsigned int verbosity=0xff;
double mineps,epsfac;
char *column=NULL;
char eps_set=0,euclidean=0,dimset=0,resize_eps;
char *outfile=NULL,stdo=1;
char *infile=NULL;

double *d_min,*d_max,d_max_max;
double **series,**delta,**corr;
double *metric,trace;
long **box,*list;
unsigned long *flist;
int emb_offset;
unsigned int ibox=BOX-1;
unsigned int *index_comp,*index_embed;

/*these are global to save time*/
int *sorted;
double *av,**mat,*matarray,*eig;

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
  fprintf(stderr,"\t-c column to read [Default: 1,..,# of components]\n");
  fprintf(stderr,"\t-m # of components,embedding dimension [Default: 1,5]\n");
  fprintf(stderr,"\t-d delay [Default: 1]\n");
  fprintf(stderr,"\t-q dimension to project to [Default: 2]\n");
  fprintf(stderr,"\t-k minimal number of neighbours [Default: 50]\n");
  fprintf(stderr,"\t-r minimal neighbourhood size \n\t\t"
	  "[Default: (interval of data)/1000]\n");
  fprintf(stderr,"\t-i # of iterations [Default: 1]\n");
  fprintf(stderr,"\t-2 use euklidean metric [Default: non euklidean]\n");
  fprintf(stderr,"\t-o name of output file \n\t\t"
	  "[Default: 'datafile'.opt.n, where n is the iteration.\n\t\t"
	  " If no -o is given, the last iteration is also"
	  " written to stdout]\n");
  fprintf(stderr,"\t-V verbosity level [Default: 7]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n\t\t"
          "2='+ average correction and trend'\n\t\t"
	  "4='+ how many points for which epsilon'\n");
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
  if ((out=check_option(in,n,'c','s')) != NULL) {
    column=out;
    dimset=1;
  }
  if ((out=check_option(in,n,'m','2')) != NULL)
    sscanf(out,"%u,%u",&comp,&embed);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(in,n,'q','u')) != NULL)
    sscanf(out,"%u",&qdim);
  if ((out=check_option(in,n,'k','u')) != NULL)
    sscanf(out,"%u",&minn);
  if ((out=check_option(in,n,'r','f')) != NULL) {
    eps_set=1;
    sscanf(out,"%lf",&mineps);
  }
  if ((out=check_option(in,n,'i','u')) != NULL)
    sscanf(out,"%u",&iterations);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'2','n')) != NULL)
    euclidean=1;
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void sort(double *x,int *n)
{
  long i,j,iswap;
  double dswap;
  
  for (i=0;i<dim;i++)
    n[i]=i;
  
  for (i=0;i<dim-1;i++)
    for (j=i+1;j<dim;j++)
      if (x[j] > x[i]) {
	dswap=x[i];
	x[i]=x[j];
	x[j]=dswap;
	iswap=n[i];
	n[i]=n[j];
	n[j]=iswap;
      }
}

void mmb(double eps)
{
  long i,x,y;
  double ieps=1.0/eps;

  for (x=0;x<BOX;x++)
    for (y=0;y<BOX;y++)
      box[x][y] = -1;
  
  for (i=emb_offset;i<length;i++) {
    x=(int)(series[0][i]*ieps)&ibox;
    y=(int)(series[comp-1][i-emb_offset]*ieps)&ibox;
    list[i]=box[x][y];
    box[x][y]=i;
  }
}

unsigned long fmn(long which,double eps)
{
  unsigned long nf=0;
  long i,i1,i2,j,j1,k,k1,li;
  long element;
  double dx=0.0;
  
  i=(int)(series[0][which]/eps)&ibox;
  j=(int)(series[comp-1][which-emb_offset]/eps)&ibox;
  
  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1&ibox;
    for (j1=j-1;j1<=j+1;j1++) {
      element=box[i2][j1&ibox];
      while (element != -1) {
	for (k=0;k<embed;k++) {
	  k1= -k*(int)delay;
	  for (li=0;li<comp;li++) {
	    dx=fabs(series[li][which+k1]-series[li][element+k1]);
	    if (dx > eps)
	      break;
	  }
	  if (dx > eps)
	    break;
	}
	if (dx <= eps)
	  flist[nf++]=element;
	element=list[element];
      }
    }
  }
  return nf;
}

void make_correction(unsigned long n,unsigned long nf)
{
  long i,i1,i2,j,j1,j2,k,k1,k2,hs;
  double help;
  
  for (i=0;i<dim;i++) {
    i1=index_comp[i];
    i2=index_embed[i];
    help=0.0;
    for (j=0;j<nf;j++)
      help += series[i1][flist[j]-i2];
    av[i]=help/nf;
  }

  for (i=0;i<dim;i++) {
    i1=index_comp[i];
    i2=index_embed[i];
    for (j=i;j<dim;j++) {
      help=0.0;
      j1=index_comp[j];
      j2=index_embed[j];
      for (k=0;k<nf;k++) {
	hs=flist[k];
	help += series[i1][hs-i2]*series[j1][hs-j2];
      }
      mat[i][j]=(help/nf-av[i]*av[j])*metric[i]*metric[j];
      mat[j][i]=mat[i][j];
    }
  }

  eigen(mat,(unsigned long)dim,eig);
  sort(eig,sorted);

  for (i=0;i<dim;i++) {
    help=0.0;
    for (j=qdim;j<dim;j++) {
      hs=sorted[j];
      for (k=0;k<dim;k++) {
	k1=index_comp[k];
	k2=index_embed[k];
	help += (series[k1][n-k2]-av[k])*mat[k][hs]*mat[i][hs]*metric[k];
      }
    }
    corr[n][i]=help/metric[i];
  }
}

void handle_trend(unsigned long n,unsigned long nf)
{
  long i,i1,i2,j;
  double help;
  
  for (i=0;i<dim;i++) {
    help=0.0;
    for (j=0;j<nf;j++)
      help += corr[flist[j]][i];
    av[i]=help/nf;
  }

  for (i=0;i<dim;i++) {
    i1=index_comp[i];
    i2=index_embed[i];
    delta[i1][n-i2] += (corr[n][i]-av[i])/(trace*metric[i]);
  }
}

void set_correction(void)
{
  long i,j;
  double *hav,*hsigma,help;

  check_alloc(hav=(double*)malloc(sizeof(double)*comp));
  check_alloc(hsigma=(double*)malloc(sizeof(double)*comp));
  for (j=0;j<comp;j++)
    hav[j]=hsigma[j]=0.0;

  for (i=0;i<length;i++)
    for (j=0;j<comp;j++) {
      hav[j] += (help=delta[j][i]);
      hsigma[j] += help*help;
    }

  for (j=0;j<comp;j++) {
    hav[j] /= length;
    hsigma[j]=sqrt(fabs(hsigma[j]/length-hav[j]*hav[j]));
  }
  if (verbosity&(VER_USR1|VER_USR2)) {
    for (i=0;i<comp;i++) {
      fprintf(stderr,"Average shift of component %ld = %e\n",i+1,
	      hav[i]*d_max[i]);
      fprintf(stderr,"Average rms correction of comp. %ld = %e\n\n",
	      i+1,hsigma[i]*d_max[i]);
    }
  }
  for (i=0;i<length;i++)
    for (j=0;j<comp;j++)
      series[j][i] -= delta[j][i];

  if (resize_eps) {
    mineps /= epsfac;
    if (verbosity&VER_USR2)
      fprintf(stderr,"Reset minimal neighbourhood size to %e\n",
	      mineps*d_max_max);
  }

  resize_eps=0;
  free(hav);
  free(hsigma);
}

int main(int argc,char **argv)
{
  char stdi=0;
  int iter,epscount,*ok;
  long i,j;
  char all_done;
  char *ofname;
  unsigned long nfound,n,allfound;
  double epsilon;
  double **hser;
  FILE *file;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  dim=comp*embed;
  emb_offset=(embed-1)*delay;

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      check_alloc(ofname=(char*)calloc(strlen(infile)+9,(size_t)1));
      sprintf(outfile,"%s.opt",infile);
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      check_alloc(ofname=(char*)calloc((size_t)14,(size_t)1));
      sprintf(outfile,"stdin.opt");
    }
  }
  else
    check_alloc(ofname=(char*)calloc(strlen(outfile)+10,(size_t)1));
  
  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&comp,"",
				     dimset,verbosity);
  else 
    series=(double**)get_multi_series(infile,&length,exclude,&comp,column,
				      dimset,verbosity);

  if (length < minn) {
    fprintf(stderr,"With %lu data you will never find %u neighbors."
	    " Exiting!\n",length,minn);
    exit(GHKSS__TOO_MANY_NEIGHBORS);
  }

  check_alloc(d_min=(double*)malloc(sizeof(double)*comp));
  check_alloc(d_max=(double*)malloc(sizeof(double)*comp));
  d_max_max=0.0;
  for (i=0;i<comp;i++) {
    rescale_data(series[i],length,&d_min[i],&d_max[i]);
    if (d_max[i] > d_max_max)
      d_max_max=d_max[i];
  }

  if (!eps_set)
    mineps=1./1000.;
  else
    mineps /= d_max_max;
  epsfac=sqrt(2.0);

  check_alloc(box=(long**)malloc(sizeof(long*)*BOX));
  for (i=0;i<BOX;i++)
    check_alloc(box[i]=(long*)malloc(sizeof(long)*BOX));

  check_alloc(list=(long*)malloc(sizeof(long)*length));
  check_alloc(flist=(unsigned long*)malloc(sizeof(long)*length));

  check_alloc(metric=(double*)malloc(sizeof(double)*dim));
  trace=0.0;
  if (euclidean) {
    for (i=0;i<dim;i++) {
      metric[i]=1.0;
      trace += 1./metric[i];
    }
  }
  else {
    for (i=0;i<dim;i++) {
      if ((i >= comp) && (i < ((long)dim-(long)comp))) 
	metric[i]=1.0;
      else 
	metric[i]=1.0e3;
      trace += 1./metric[i];
    }
  }

  check_alloc(corr=(double**)malloc(sizeof(double*)*length));
  for (i=0;i<length;i++)
    check_alloc(corr[i]=(double*)malloc(sizeof(double)*dim));
  check_alloc(ok=(int*)malloc(sizeof(int)*length));
  check_alloc(delta=(double**)malloc(sizeof(double*)*comp));
  for (i=0;i<comp;i++)
    check_alloc(delta[i]=(double*)malloc(sizeof(double)*length));
  check_alloc(index_comp=(unsigned int*)malloc(sizeof(int)*dim));
  check_alloc(index_embed=(unsigned int*)malloc(sizeof(int)*dim));
  check_alloc(av=(double*)malloc(sizeof(double)*dim));
  check_alloc(sorted=(int*)malloc(sizeof(int)*dim));
  check_alloc(eig=(double*)malloc(sizeof(double)*dim));
  check_alloc(matarray=(double*)malloc(sizeof(double)*dim*dim));
  check_alloc(mat=(double**)malloc(sizeof(double*)*dim));
  for (i=0;i<dim;i++)
    mat[i]=(double*)(matarray+dim*i);
  check_alloc(hser=(double**)malloc(sizeof(double*)*comp));

  for (i=0;i<dim;i++) {
    index_comp[i]=i%comp;
    index_embed[i]=(i/comp)*delay;
  }

  resize_eps=0;
  for (iter=1;iter<=iterations;iter++) {
    for (i=0;i<length;i++) {
      ok[i]=0;
      for (j=0;j<dim;j++)
	corr[i][j]=0.0;
      for (j=0;j<comp;j++)
	delta[j][i]=0.0;
    }
    epsilon=mineps;
    all_done=0;
    epscount=1;
    allfound=0;
    if (verbosity&(VER_USR1|VER_USR2))
      fprintf(stderr,"Starting iteration %d\n",iter);
    while(!all_done) {
      mmb(epsilon);
      all_done=1;
      for (n=emb_offset;n<length;n++)
	if (!ok[n]) {
	  nfound=fmn(n,epsilon);
	  if (nfound >= minn) {
	    make_correction(n,nfound);
	    ok[n]=epscount;
	    if (epscount == 1)
	      resize_eps=1;
	    allfound++;
	  }
	  else
	    all_done=0;
	}
      if (verbosity&VER_USR2)
	fprintf(stderr,"Corrected %ld points with epsilon= %e\n",allfound,
		epsilon*d_max_max);
      epsilon *= epsfac;
      epscount++;
    }
    if (verbosity&VER_USR2)
      fprintf(stderr,"Start evaluating the trend\n");

    epsilon=mineps;
    allfound=0;
    for (i=1;i<epscount;i++) {
      mmb(epsilon);
      for (n=emb_offset;n<length;n++)
	if (ok[n] == i) {
	  nfound=fmn(n,epsilon);
	  handle_trend(n,nfound);
	  allfound++;
	}
      if (verbosity&VER_USR2)
	fprintf(stderr,"Trend subtracted for %ld points with epsilon= %e\n",
		allfound,epsilon*d_max_max);
      epsilon *= epsfac;
    }
    set_correction();
    
    sprintf(ofname,"%s.%d",outfile,iter);
    test_outfile(ofname);

    file=fopen(ofname,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n\n",ofname);
    for (i=0;i<length;i++) {
      for (j=0;j<comp;j++) {
	fprintf(file,"%e ",series[j][i]*d_max[j]+d_min[j]);
      }
      fprintf(file,"\n");
      if (stdo && (iter == iterations)) {
	for (j=0;j<comp;j++)
	  fprintf(stdout,"%e ",series[j][i]*d_max[j]+d_min[j]);
	fprintf(stdout,"\n");
      }
    }
    fclose(file);
  }

  return 0;
}
