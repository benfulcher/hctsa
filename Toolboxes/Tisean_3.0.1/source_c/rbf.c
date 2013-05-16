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
/*Author: Rainer Hegger. Last modified: Mar 11, 2002 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"
#include <math.h>

#define WID_STR "Fits a RBF-model to the data"

char *outfile=NULL,stdo=1,MAKECAST=0;
char *infile=NULL;
char setdrift=1;
int DIM=2,DELAY=1,CENTER=10,STEP=1;
unsigned int COLUMN=1;
unsigned int verbosity=0xff;
long CLENGTH=1000;
unsigned long LENGTH=ULONG_MAX,INSAMPLE=ULONG_MAX,exclude=0;

double *series,*coefs;
double varianz,interval,min;
double **center;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin\n");
  fprintf(stderr,"\t-l # of data to use [default: all from file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default: 0]\n");
  fprintf(stderr,"\t-c column to read [default: %u]\n",COLUMN);
  fprintf(stderr,"\t-m embedding dimension [default: %d]\n",DIM);
  fprintf(stderr,"\t-d delay [default: %d]\n",DELAY);
  fprintf(stderr,"\t-p number of centers [default: %d]\n",CENTER);
  fprintf(stderr,"\t-X deactivate drift [default: activated]\n");
  fprintf(stderr,"\t-s steps to forecast [default: %d]\n",STEP);
  fprintf(stderr,"\t-n # of points for insample [default: # of data]\n");
  fprintf(stderr,"\t-L steps to cast [default: none]\n");
  fprintf(stderr,"\t-o output file name [default: 'datafile'.rbf]\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
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
  if ((out=check_option(in,n,'p','u')) != NULL)
    sscanf(out,"%u",&CENTER);
  if ((out=check_option(in,n,'X','n')) != NULL)
    setdrift=0;
  if ((out=check_option(in,n,'s','u')) != NULL)
    sscanf(out,"%u",&STEP);
  if ((out=check_option(in,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(in,n,'n','u')) != NULL)
    sscanf(out,"%lu",&INSAMPLE);
  if ((out=check_option(in,n,'L','u')) != NULL) {
    MAKECAST=1;
    sscanf(out,"%lu",&CLENGTH);
  }
  if ((out=check_option(in,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

double avdistance(void)
{
  int i,j,k;
  double dist=0.0;
  
  for (i=0;i<CENTER;i++)
    for (j=0;j<CENTER;j++)
      if (i != j)
	for (k=0;k<DIM;k++)
	  dist += sqr(center[i][k]-center[j][k]);

  return sqrt(dist/(CENTER-1)/CENTER/DIM);
}

double rbf(double *act,double *cen)
{
  static double denum;
  double r=0;
  int i;

  denum=2.0*varianz*varianz;

  for (i=0;i<DIM;i++)
    r += sqr(*(act-i*DELAY)-cen[i]);
  
  return exp(-r/denum);
}

void drift(void) 
{
  double *force,h,h1,step=1e-2,step1;
  int i,j,k,l,d2=DIM;

  check_alloc(force=(double*)malloc(sizeof(double)*d2));
  for (l=0;l<20;l++) {
    for (i=0;i<CENTER;i++) {
      for (j=0;j<d2;j++) {
        force[j]=0.0;
        for (k=0;k<CENTER;k++) {
          if (k != i) {
            h=center[i][j]-center[k][j];
            force[j] += h/sqr(h)/fabs(h);
          }
        }
      }
      h=0.0;
      for (j=0;j<d2;j++) 
        h += sqr(force[j]);
      step1=step/sqrt(h);
      for (j=0;j<d2;j++) {
        h1 = step1*force[j];
        if (((center[i][j]+h1) > -0.1) && ((center[i][j]+h1) < 1.1))
          center[i][j] += h1;
      }
    }
  }
  free(force);
}

void make_fit(void)
{
  double **mat,*hcen;
  double h;
  int i,j,n,nst;

  check_alloc(mat=(double**)malloc(sizeof(double*)*(CENTER+1)));
  for (i=0;i<=CENTER;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*(CENTER+1)));
  check_alloc(hcen=(double*)malloc(sizeof(double)*CENTER));

  for (i=0;i<=CENTER;i++) {
    coefs[i]=0.0;
    for (j=0;j<=CENTER;j++)
      mat[i][j]=0.0;
  }

  for (n=(DIM-1)*DELAY;n<INSAMPLE-STEP;n++) {
    nst=n+STEP;
    for (i=0;i<CENTER;i++)
      hcen[i]=rbf(&series[n],center[i]);
    coefs[0] += series[nst];
    mat[0][0] += 1.0;
    for (i=1;i<=CENTER;i++)
      mat[i][0] += hcen[i-1];
    for (i=1;i<=CENTER;i++) {
      coefs[i] += series[nst]*(h=hcen[i-1]);
      for (j=1;j<=i;j++)
	mat[i][j] += h*hcen[j-1];
    }
  }
  
  h=(double)(INSAMPLE-STEP-(DIM-1)*DELAY);
  for (i=0;i<=CENTER;i++) {
    coefs[i] /= h;
    for (j=0;j<=i;j++) {
      mat[i][j] /= h;
      mat[j][i]=mat[i][j];
    }
  }

  solvele(mat,coefs,(unsigned int)(CENTER+1));

  for (i=0;i<=CENTER;i++)
    free(mat[i]);
  free(mat);
  free(hcen);
}

double forecast_error(unsigned long i0,unsigned long i1)
{
  int i,n;
  double h,error=0.0;

  for (n=i0+(DIM-1)*DELAY;n<i1-STEP;n++) {
    h=coefs[0];
    for (i=1;i<=CENTER;i++)
      h += coefs[i]*rbf(&series[n],center[i-1]);
    error += (series[n+STEP]-h)*(series[n+STEP]-h);
  }
  
  return sqrt(error/(i1-i0-STEP-(DIM-1)*DELAY));
}

void make_cast(FILE *out)
{
  double *cast,new_el;
  int i,n,dim;
  
  dim=(DIM-1)*DELAY;
  check_alloc(cast=(double*)malloc(sizeof(double)*(dim+1)));
  for (i=0;i<=dim;i++)
    cast[i]=series[LENGTH-1-dim+i];
  
  for (n=0;n<CLENGTH;n++) {
    new_el=coefs[0];
    for (i=1;i<=CENTER;i++)
      new_el += coefs[i]*rbf(&cast[dim],center[i-1]);
    fprintf(out,"%e\n",new_el*interval+min);
    for (i=0;i<dim;i++)
      cast[i]=cast[i+1];
    cast[dim]=new_el;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  int i,j,cstep;
  double sigma,av;
  FILE *file=NULL;

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
      strcat(outfile,".rbf");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.rbf");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  series=(double*)get_series(infile,&LENGTH,exclude,COLUMN,verbosity);
  rescale_data(series,LENGTH,&min,&interval);
  variance(series,LENGTH,&av,&varianz);

  if (INSAMPLE > LENGTH)
    INSAMPLE=LENGTH;
  
  if (CENTER > LENGTH) 
    CENTER = LENGTH;
  
  if (MAKECAST)
    STEP=1;
  
  check_alloc(coefs=(double*)malloc(sizeof(double)*(CENTER+1)));
  check_alloc(center=(double**)malloc(sizeof(double*)*CENTER));
  for (i=0;i<CENTER;i++)
    check_alloc(center[i]=(double*)malloc(sizeof(double)*DIM));
  
  cstep=LENGTH-1-(DIM-1)*DELAY;
  for (i=0;i<CENTER;i++)
    for (j=0;j<DIM;j++)
      center[i][j]=series[(DIM-1)*DELAY-j*DELAY+(i*cstep)/(CENTER-1)];

  if (setdrift)
    drift();
  varianz=avdistance();
  make_fit();

  if (!stdo) {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(file,"#Center points used:\n");
    for (i=0;i<CENTER;i++) {
      fprintf(file,"#");
      for (j=0;j<DIM;j++)
	fprintf(file," %e",center[i][j]*interval+min);
      fprintf(file,"\n");
    }
    fprintf(file,"#variance= %e\n",varianz*interval);
    fprintf(file,"#Coefficients:\n");
    fprintf(file,"#%e\n",coefs[0]*interval+min);
    for (i=1;i<=CENTER;i++)
      fprintf(file,"#%e\n",coefs[i]*interval);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"#Center points used:\n");
    for (i=0;i<CENTER;i++) {
      fprintf(stdout,"#");
      for (j=0;j<DIM;j++)
	fprintf(stdout," %e",center[i][j]*interval+min);
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"#variance= %e\n",varianz*interval);
    fprintf(stdout,"#Coefficients:\n");
    fprintf(stdout,"#%e\n",coefs[0]*interval+min);
    for (i=1;i<=CENTER;i++)
      fprintf(stdout,"#%e\n",coefs[i]*interval);
  }
  av=sigma=0.0;
  for (i=0;i<INSAMPLE;i++) {
    av += series[i];
    sigma += series[i]*series[i];
  }
  av /= INSAMPLE;
  sigma=sqrt(fabs(sigma/INSAMPLE-av*av));
  if (!stdo)
    fprintf(file,"#insample error= %e\n",forecast_error(0LU,INSAMPLE)/sigma);
  else
    fprintf(stdout,"#insample error= %e\n",forecast_error(0LU,INSAMPLE)/sigma);

  if (INSAMPLE < LENGTH) {
    av=sigma=0.0;
    for (i=INSAMPLE;i<LENGTH;i++) {
      av += series[i];
      sigma += series[i]*series[i];
    }
    av /= (LENGTH-INSAMPLE);
    sigma=sqrt(fabs(sigma/(LENGTH-INSAMPLE)-av*av));
    if (!stdout)
      fprintf(file,"#out of sample error= %e\n",
	      forecast_error(INSAMPLE,LENGTH)/sigma);
    else
      fprintf(stdout,"#out of sample error= %e\n",
	      forecast_error(INSAMPLE,LENGTH)/sigma);
  }

  if (MAKECAST) {
    if (!stdo)
      make_cast(file);
    else
      make_cast(stdout);
  }

  if (!stdo)
    fclose(file);

  return 0;
}
