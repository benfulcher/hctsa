/*Author: Rainer Hegger Last modified: Sep 4, 1999 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "routines/tsa.h"

#define WID_STR "Performs a global SVD"

unsigned long LENGTH=ULONG_MAX,exclude=0;
unsigned int DIM=2,LDIM=2,COLUMN=1,DELAY=1;
unsigned int verbosity=0xff;
char *outfile=NULL,stout=1,make_base=0,dim_set=0;
char *infile=NULL;
double *series,av=0.0;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to use [Default: whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignore [Default: 0]\n");
  fprintf(stderr,"\t-c column to read [Default: 1]\n");
  fprintf(stderr,"\t-m dimension to use [Default: 2]\n");
  fprintf(stderr,"\t-d delay to use [Default: 1]\n");
  fprintf(stderr,"\t-q write also the vectors in the SVD basis after\n"
	  "\t\tprojecting onto # dimensions [Default: none]\n");
  fprintf(stderr,"\t-o output file name \n\t\t[Default: stdout; -o without "
	  "value means 'datafile'.svd]\n");
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
    sscanf(out,"%lu",&LENGTH);
  if ((out=check_option(in,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(in,n,'c','u')) != NULL)
    sscanf(out,"%u",&COLUMN);
  if ((out=check_option(in,n,'m','u')) != NULL)
    sscanf(out,"%u",&DIM);
  if ((out=check_option(in,n,'d','u')) != NULL)
    sscanf(out,"%u",&DELAY);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'q','u')) != NULL) {
    make_base=1;
    sscanf(out,"%u",&LDIM);
    if (LDIM < 1) LDIM=1;
    if (LDIM > DIM) LDIM=DIM;
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void ordne(double *lyap,int *ord)
{
  long i,j,maxi;
  double max;
  
  for (i=0;i<DIM;i++)
    ord[i]=i;

  for (i=0;i<DIM-1;i++)
    for (j=i+1;j<DIM;j++)
      if (lyap[i] < lyap[j]) {
	max=lyap[i];
	lyap[i]=lyap[j];
	lyap[j]=max;
	maxi=ord[i];
	ord[i]=ord[j];
	ord[j]=maxi;
      }
}

void make_pca(void)
{
  unsigned int i,j,k,l,hi,hj;
  int *ord;
  double **mat,*eig,*off,*vec,help=0.0;
  FILE *fout=NULL;

  check_alloc(ord=(int*)malloc(sizeof(int)*DIM));
  check_alloc(off=(double*)malloc(sizeof(double)*DIM));
  check_alloc(eig=(double*)malloc(sizeof(double)*DIM));
  check_alloc(vec=(double*)malloc(sizeof(double)*DIM));
  check_alloc(mat=(double**)malloc(sizeof(double*)*DIM));
  for (i=0;i<DIM;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*DIM));

  
  for (i=0;i<DIM;i++) {
    hi=i*DELAY;
    for (j=i;j<DIM;j++) {
      hj=j*DELAY;
      mat[i][j]=0.0;
      for (k=0;k<LENGTH-(DIM-1)*DELAY;k++)
	mat[i][j] += series[k+hi]*series[k+hj];
      mat[j][i]=(mat[i][j] /= (LENGTH-DIM));
    }
  }

  eig1(mat,(long)DIM,eig,off);
  eig2(eig,off,(long)DIM,mat);
  ordne(eig,ord);
  
  if (!stout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
  }

  for (i=0;i<DIM;i++)
    if (!make_base) {
      if (stout)
	fprintf(stdout,"%d %e\n",i,eig[i]);
      else
	fprintf(fout,"%d %e\n",i,eig[i]);
    }
    else {
      if (stout)
	fprintf(stdout,"#%d %e\n",i,eig[i]);
      else
	fprintf(fout,"#%d %e\n",i,eig[i]);
    }

  
  if (make_base) {
    if (LDIM == DIM) {
      for (i=0;i<=LENGTH-(DIM-1)*DELAY;i++) {
	for (j=0;j<DIM;j++)
	  vec[j]=0.0;
	for (j=0;j<LDIM;j++) {
	  hj=(unsigned int)ord[j];
	  for (k=0;k<DIM;k++)
	    vec[j] += series[i+k*DELAY]*mat[k][hj];
	}
	for (j=0;j<LDIM;j++)
	  if (stout)
	    fprintf(stdout,"%e ",vec[j]+av);
	  else
	    fprintf(fout,"%e ",vec[j]+av);
	if (stout)
	  fprintf(stdout,"\n");
	else
	  fprintf(fout,"\n");
      }
    }
    else {
      for (i=0;i<=LENGTH-(DIM-1)*DELAY;i++) {
	for (l=0;l<DIM;l++) {
	  help=0.0;
	  for (j=0;j<LDIM;j++) {
	    hj=(unsigned int)ord[j];
	    for (k=0;k<DIM;k++)
	      help += series[i+k*DELAY]*mat[k][hj]*mat[l][hj];
	  }
	}
	if (stout)
	  fprintf(stdout,"%e\n",help+av);
	else
	  fprintf(fout,"%e\n",help+av);
      }
    }
  }
  if (!stout)
    fclose(fout);
}

int main(int argc,char **argv)
{
  char stdi=0;
  int i;
  double rms;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,&COLUMN,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".svd");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.svd");
    }
  }
  if (!stout)
    test_outfile(outfile);

  series=(double*)get_series(infile,&LENGTH,exclude,COLUMN,verbosity);
  variance(series,LENGTH,&av,&rms);
  for (i=0;i<LENGTH;i++)
    series[i] -= av;

  make_pca();

  return 0;
}
