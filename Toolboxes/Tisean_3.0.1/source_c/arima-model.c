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
/*Author: Rainer Hegger, Last modified: Feb 6, 2006 */
/*Changes:
  Feb 4, 2006: First version
  Feb 6, 2006: Find and remove bugs (1)
  Feb 11, 2006: Add rand_arb_dist to iterate_***_model
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "routines/tsa.h"

#define WID_STR "Fits an multivariate ARIMA model to the data and gives\
 the coefficients\n\tand the residues (or an iterated model)"

unsigned long length=ULONG_MAX,exclude=0;
unsigned int dim=1,poles=10,ilength,ITER=50;
unsigned int arpoles=0,ipoles=0,mapoles=0,offset;
unsigned int verbosity=1;
char *outfile=NULL,*column=NULL,stdo=1,dimset=0,run_model=0,arimaset=0;
char *infile=NULL;
double **series,convergence=1.0e-3;

double *my_average;
unsigned long ardim,armadim;
unsigned int **aindex;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
	  " as a possible"
	  " datafile.\nIf no datafile is given stdin is read. Just - also"
	  " means stdin\n");
  fprintf(stderr,"\t-l length of file [default is whole file]\n");
  fprintf(stderr,"\t-x # of lines to be ignored [default is 0]\n");
  fprintf(stderr,"\t-m dimension [default is 1]\n");
  fprintf(stderr,"\t-c columns to read [default is 1,...,dimension]\n");
  fprintf(stderr,"\t-p order of initial AR-Fit [default is %u]\n",poles);
  fprintf(stderr,"\t-P order of AR,I,MA-Fit [default is %u,%u,%u]\n",
	  arpoles,ipoles,mapoles);
  fprintf(stderr,"\t-I # of arima iterations [default is %u]\n",ITER);
  fprintf(stderr,"\t-e accuracy of convergence [default is %lf]\n",convergence);
  fprintf(stderr,"\t-s length of iterated model [default no iteration]\n");
  fprintf(stderr,"\t-o output file name [default is 'datafile'.ari]\n");
  fprintf(stderr,"\t-V verbosity level [default is 1]\n\t\t"
	  "0='only panic messages'\n\t\t"
	  "1='+ input/output messages'\n\t\t"
	  "2='+ print residuals though iterating a model'\n\t\t"
	  "4='+ print original data plus residuals'\n");
  fprintf(stderr,"\t-h show these options\n\n");
  exit(0);
}

void scan_options(int argc,char **argv)
{
  char *out;

  if ((out=check_option(argv,argc,'p','u')) != NULL) {
    sscanf(out,"%u",&poles);
    if (poles < 1) {
      fprintf(stderr,"The order should at least be one!\n");
      exit(127);
    }
  }
  if ((out=check_option(argv,argc,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(argv,argc,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(argv,argc,'m','u')) != NULL) {
    sscanf(out,"%u",&dim);
    dimset=1;
  }
  if ((out=check_option(argv,argc,'P','3')) != NULL) {
    sscanf(out,"%u,%u,%u",&arpoles,&ipoles,&mapoles);
    if ((arpoles+ipoles+mapoles)>0)
      arimaset=1;
  }
  if ((out=check_option(argv,argc,'I','u')) != NULL)
    sscanf(out,"%u",&ITER);
  if ((out=check_option(argv,argc,'e','f')) != NULL)
    sscanf(out,"%lf",&convergence);
  if ((out=check_option(argv,argc,'c','u')) != NULL)
    column=out;
  if ((out=check_option(argv,argc,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,argc,'s','u')) != NULL) {
    sscanf(out,"%u",&ilength);
    run_model=1;
  }
  if ((out=check_option(argv,argc,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void make_difference(void)
{
  unsigned long i,d;

  for (i=length-1;i>0;i--)
    for (d=0;d<dim;d++)
      series[d][i]=series[d][i]-series[d][i-1];
}

unsigned int** make_ar_index(void)
{
  unsigned int** ar_index;
  unsigned long i;

  check_alloc(ar_index=(unsigned int**)malloc(sizeof(unsigned int*)*2));
  for (i=0;i<2;i++)
    check_alloc(ar_index[i]=(unsigned int*)
		malloc(sizeof(unsigned int)*ardim));
  for (i=0;i<ardim;i++) {
    ar_index[0][i]=i/poles;
    ar_index[1][i]=i%poles;
  }
  return ar_index;
}

unsigned int** make_arima_index(unsigned int ars,unsigned int mas)
{
  unsigned int** arima_index;
  unsigned int armad;
  unsigned long i,i0;

  armad=(ars+mas)*dim;
  check_alloc(arima_index=(unsigned int**)malloc(sizeof(unsigned int*)*2));
  for (i=0;i<2;i++)
    check_alloc(arima_index[i]=(unsigned int*)
		malloc(sizeof(unsigned int)*armad));
  for (i=0;i<ars*dim;i++) {
    arima_index[0][i]=i/ars;
    arima_index[1][i]=i%ars;
  }
  i0=ars*dim;
  for (i=0;i<mas*dim;i++) {
    arima_index[0][i+i0]=dim+i/mas;
    arima_index[1][i+i0]=i%mas;
  }

  return arima_index;
}

void set_averages_to_zero(void)
{
  double var;
  long i,j;
  
  for (i=0;i<dim;i++) {
    variance(series[i],length,&my_average[i],&var);
    for (j=0;j<length;j++)
      series[i][j] -= my_average[i];
  }
}

double** build_matrix(double **mat,unsigned int size)
{
  long n,i,j,is,id,js,jd;
  double norm;
  
  norm=1./((double)length-1.0-(double)poles-(double)offset);

  for (i=0;i<size;i++) {
    id=aindex[0][i];
    is=aindex[1][i];
    for (j=i;j<size;j++) {
      jd=aindex[0][j];
      js=aindex[1][j];
      mat[i][j]=0.0;
      for (n=offset+poles-1;n<length-1;n++)
	mat[i][j] += series[id][n-is]*series[jd][n-js];
      mat[i][j] *= norm;
      mat[j][i]=mat[i][j];
    }
  }

  return invert_matrix(mat,size);
}

void build_vector(double *vec,unsigned int size,long comp)
{
  long i,is,id,n;
  double norm;

  norm=1./((double)length-1.0-(double)poles-(double)offset);

  for (i=0;i<size;i++) {
    id=aindex[0][i];
    is=aindex[1][i];
    vec[i]=0.0;
    for (n=offset+poles-1;n<length-1;n++)
      vec[i] += series[comp][n+1]*series[id][n-is];
    vec[i] *= norm;
  }
}

double* multiply_matrix_vector(double **mat,double *vec,unsigned int size)
{
  long i,j;
  double *new_vec;

  check_alloc(new_vec=(double*)malloc(sizeof(double)*size));

  for (i=0;i<size;i++) {
    new_vec[i]=0.0;
    for (j=0;j<size;j++)
      new_vec[i] += mat[i][j]*vec[j];
  }

  return new_vec;
}

double* make_residuals(double **diff,double **coeff,unsigned int size)
{
  long n,n1,d,i,is,id;
  double *resi;
  
  check_alloc(resi=(double*)malloc(sizeof(double)*dim));
  for (i=0;i<dim;i++)
    resi[i]=0.0;

  for (n=poles-1;n<length-1;n++) {
    n1=n+1;
    for (d=0;d<dim;d++) {
      diff[d][n1]=series[d][n1];
      for (i=0;i<size;i++) {
	id=aindex[0][i];
	is=aindex[1][i];
	diff[d][n1] -= coeff[d][i]*series[id][n-is];
      }
      resi[d] += sqr(diff[d][n1]);
    }
  }

  for (i=0;i<dim;i++)
    resi[i]=sqrt(resi[i]/((double)length-(double)poles));

  return resi;
}

void iterate_model(double **coeff,double *sigma,double **diff,FILE *file)
{
  long i,j,i1,i2,n,d;
  double **iterate,*swap,**myrand;
  
  check_alloc(iterate=(double**)malloc(sizeof(double*)*(poles+1)));
  for (i=0;i<=poles;i++)
    check_alloc(iterate[i]=(double*)malloc(sizeof(double)*dim));

  check_alloc(myrand=(double**)malloc(sizeof(double*)*dim));
  for (i=0;i<dim;i++)
    myrand[i]=rand_arb_dist(diff[i],length,ilength+poles,100,0x44325);

  rnd_init(0x44325);
  for (i=0;i<1000;i++)
    rnd_long();
  for (i=0;i<dim;i++)
    for (j=0;j<poles;j++)
      iterate[j][i]=myrand[i][j];
  
  for (n=0;n<ilength;n++) {
    for (d=0;d<dim;d++) {
      iterate[poles][d]=myrand[d][n+poles];
      for (i1=0;i1<dim;i1++)
	for (i2=0;i2<poles;i2++)
	  iterate[poles][d] += coeff[d][i1*poles+i2]*iterate[poles-1-i2][i1];
    }
    if (file != NULL) {
      for (d=0;d<dim;d++)
	fprintf(file,"%e ",iterate[poles][d]);
      fprintf(file,"\n");
    }
    else {
      for (d=0;d<dim;d++)
	printf("%e ",iterate[poles][d]);
      printf("\n");
    }

    swap=iterate[0];
    for (i=0;i<poles;i++)
      iterate[i]=iterate[i+1];
    iterate[poles]=swap;
  }

  for (i=0;i<=poles;i++)
    free(iterate[i]);
  free(iterate);

  for (i=0;i<dim;i++)
    free(myrand[i]);
  free(myrand);
}

void iterate_arima_model(double **coeff,double *sigma,double **diff,FILE *file)
{
  double **iterate,*swap,**myrand;
  unsigned long i,j,n,is,id;

  check_alloc(iterate=(double**)malloc(sizeof(double*)*(poles+1)));
  for (i=0;i<=poles;i++)
    check_alloc(iterate[i]=(double*)malloc(sizeof(double)*2*dim));

  check_alloc(myrand=(double**)malloc(sizeof(double*)*dim));
  for (i=0;i<dim;i++)
    myrand[i]=rand_arb_dist(diff[i],length,ilength+poles,100,0x44325);

  rnd_init(0x44325);
  for (i=0;i<1000;i++)
    rnd_long();
  for (i=0;i<dim;i++)
    for (j=0;j<poles;j++)
      iterate[j][i]=iterate[j][dim+i]=myrand[i][j];

  for (n=0;n<ilength;n++) {
    for (i=0;i<dim;i++)
      iterate[poles][i]=iterate[poles][i+dim]=myrand[i][n+poles];

    for (j=0;j<dim;j++) {
      for (i=0;i<armadim;i++) {
	id=aindex[0][i];
	is=aindex[1][i];
	iterate[poles][j] += coeff[j][i]*iterate[poles-1-is][id];
      }
    }

    if (file != NULL) {
      for (i=0;i<dim;i++)
	fprintf(file,"%e ",iterate[poles][i]);
      fprintf(file,"\n");
    }
    else {
      for (i=0;i<dim;i++)
	printf("%e ",iterate[poles][i]);
      printf("\n");
    }

    swap=iterate[0];
    for (i=0;i<poles;i++)
      iterate[i]=iterate[i+1];
    iterate[poles]=swap;
  }

  for (i=0;i<=poles;i++)
    free(iterate[i]);
  free(iterate);
  for (i=0;i<dim;i++)
    free(myrand[i]);
  free(myrand);
}

int main(int argc,char **argv)
{
  char stdi=0;
  double *pm;
  long i,j,iter,hj,realiter=0;
  unsigned int size,is,id;
  FILE *file;
  double **mat,**inverse,*vec,**coeff,**diff,**hseries;
  double **oldcoeff,*diffcoeff=NULL;
  double hdiff,**xdiff=NULL,avpm;
  double loglikelihood,aic,alldiff;
  
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
      strcat(outfile,".ari");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.ari");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",dimset,
				      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,column,
				      dimset,verbosity);

  check_alloc(my_average=(double*)malloc(sizeof(double)*dim));

  for (i=0;i<ipoles;i++)
    make_difference();

  for (i=0;i<dim;i++)
    series[i] += ipoles;
  length -= ipoles;

  set_averages_to_zero();

  if (poles >= length) {
    fprintf(stderr,"It makes no sense to have more poles than data! Exiting\n");
    exit(AR_MODEL_TOO_MANY_POLES);
  }
  if (arimaset) {
    if ((arpoles >= length) || (mapoles >= length)) {
      fprintf(stderr,"It makes no sense to have more poles than data! Exiting\n");
      exit(AR_MODEL_TOO_MANY_POLES);
    }
  }
 
  ardim=poles*dim;
  aindex=make_ar_index();

  check_alloc(vec=(double*)malloc(sizeof(double)*ardim));
  check_alloc(mat=(double**)malloc(sizeof(double*)*ardim));
  for (i=0;i<ardim;i++)
    check_alloc(mat[i]=(double*)malloc(sizeof(double)*ardim));

  check_alloc(coeff=(double**)malloc(sizeof(double*)*dim));
  inverse=build_matrix(mat,ardim);
  for (i=0;i<dim;i++) {
    build_vector(vec,ardim,i);
    coeff[i]=multiply_matrix_vector(inverse,vec,ardim);
  }

  check_alloc(diff=(double**)malloc(sizeof(double*)*dim));
  for (i=0;i<dim;i++)
    check_alloc(diff[i]=(double*)malloc(sizeof(double)*length));

  pm=make_residuals(diff,coeff,ardim);

  free(vec);
  for (i=0;i<ardim;i++) {
    free(mat[i]);
    free(inverse[i]);
  }
  free(mat);
  free(inverse);
  size=ardim;
  
  if (arimaset) {
    offset=poles;
    for (i=0;i<2;i++)
      free(aindex[i]);
    free(aindex);

    for (i=0;i<dim;i++)
      free(coeff[i]);
    free(coeff);
    check_alloc(xdiff=(double**)malloc(sizeof(double*)*ITER));
    for (i=0;i<ITER;i++)
      check_alloc(xdiff[i]=(double*)malloc(sizeof(double)*dim));

    armadim=(arpoles+mapoles)*dim;
    aindex=make_arima_index(arpoles,mapoles);
    size=armadim;

    check_alloc(hseries=(double**)malloc(sizeof(double*)*2*dim));
    for (i=0;i<dim;i++) {
      check_alloc(hseries[i]=(double*)malloc(sizeof(double)*length));
      check_alloc(hseries[i+dim]=(double*)malloc(sizeof(double)*length));
      for (j=0;j<length;j++) {
	hseries[i][j]=series[i][j];
	hseries[i+dim][j]=diff[i][j];
      }
    }

    for (i=0;i<dim;i++)
      free(series[i]-ipoles);
    free(series);

    series=hseries;

    check_alloc(oldcoeff=(double**)malloc(sizeof(double*)*dim));
    for (i=0;i<dim;i++) {
      check_alloc(oldcoeff[i]=(double*)malloc(sizeof(double)*armadim));
      for (j=0;j<armadim;j++)
	oldcoeff[i][j]=0.0;
    }
    check_alloc(diffcoeff=(double*)malloc(sizeof(double)*ITER));

    for (iter=1;iter<=ITER;iter++) {
      check_alloc(vec=(double*)malloc(sizeof(double)*armadim));
      check_alloc(mat=(double**)malloc(sizeof(double*)*armadim));
      for (i=0;i<armadim;i++)
	check_alloc(mat[i]=(double*)malloc(sizeof(double)*armadim));

      check_alloc(coeff=(double**)malloc(sizeof(double*)*dim));

      poles=(arpoles > mapoles)? arpoles:mapoles;

      offset += poles;
      inverse=build_matrix(mat,armadim);

      for (i=0;i<dim;i++) {
	build_vector(vec,armadim,i);
	coeff[i]=multiply_matrix_vector(inverse,vec,armadim);
      }

      pm=make_residuals(diff,coeff,armadim);

      for (j=0;j<dim;j++) {
	hdiff=0.0;
	hj=j+dim;
	for (i=offset;i<length;i++)
	  hdiff += sqr(series[hj][i]-diff[j][i]);
	for (i=0;i<length;i++) {
	  series[hj][i]=diff[j][i];
	}
	xdiff[iter-1][j]=sqrt(hdiff/(double)(length-offset));
      }

      free(vec);
      for (i=0;i<armadim;i++) {
	free(mat[i]);
	free(inverse[i]);
      }
      free(mat);
      free(inverse);

      diffcoeff[iter-1]=0.0;
      for (i=0;i<dim;i++)
	for (j=0;j<dim;j++) {
	  diffcoeff[iter-1] += sqr(coeff[i][j]-oldcoeff[i][j]);
	  oldcoeff[i][j]=coeff[i][j];
	}
      diffcoeff[iter-1]=sqrt(diffcoeff[iter-1]/(double)armadim);
      alldiff=xdiff[iter-1][0];
      for (i=1;i<dim;i++)
	if (xdiff[iter-1][i] > alldiff)
	  alldiff=xdiff[iter-1][i];
      realiter=iter;
      if (alldiff < convergence)
	iter=ITER;
  
      if (iter < ITER) {
	for (i=0;i<dim;i++)
	  free(coeff[i]);
	free(coeff);
      }
    }
  }

  if (stdo) {
    if (arimaset) {
      printf("#convergence of residuals in arima fit\n");
      for (i=0;i<realiter;i++) {
	printf("#iteration %ld ",i+1);
	for (j=0;j<dim;j++)
	  printf("%e ",xdiff[i][j]);
	printf("%e",diffcoeff[i]);
	printf("\n");
      }
    }
    avpm=pm[0]*pm[0];
    loglikelihood= -log(pm[0]);
    for (i=1;i<dim;i++) {
      avpm += pm[i]*pm[i];
      loglikelihood -= log(pm[i]);
    }
    loglikelihood *= ((double)length);
    loglikelihood += -((double)length)*
      ((1.0+log(2.*M_PI))*dim)/2.0;
    avpm=sqrt(avpm/dim);
    printf("#average forcast error= %e\n",avpm);
    printf("#individual forecast errors: ");
     for (i=0;i<dim;i++)
      printf("%e ",pm[i]);
    printf("\n");
    if (arimaset)
      aic=2.0*(arpoles+mapoles)-2.0*loglikelihood;
    else
      aic=2.0*poles-2.0*loglikelihood;
    printf("#Log-Likelihood= %e\t AIC= %e\n",loglikelihood,aic);
    for (i=0;i<size;i++) {
      id=aindex[0][i];
      is=aindex[1][i];
      if (id < dim)
	printf("#x_%u(n-%u) ",id+1,is);
      else
	printf("#e_%u(n-%u) ",id+1-dim,is);
      for (j=0;j<dim;j++)
	printf("%e ",coeff[j][i]);
      printf("\n");
    }
    if (!run_model || (verbosity&VER_USR1)) {
      for (i=poles;i<length;i++) {
	if (run_model)
	  printf("#");
	for (j=0;j<dim;j++)
	  if (verbosity&VER_USR2)
	    printf("%e %e ",series[j][i]+my_average[j],diff[j][i]);
	  else
	    printf("%e ",diff[j][i]);
	printf("\n");
      }
    }
    if (run_model && (ilength > 0)) {
      if (!arimaset)
	iterate_model(coeff,pm,diff,NULL);
      else 
	iterate_arima_model(coeff,pm,diff,NULL);
    }
  }
  else {
    file=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for output\n",outfile);
    if (arimaset) {
      fprintf(file,"#convergence of residuals in arima fit\n");
      for (i=0;i<realiter;i++) {
	fprintf(file,"#iteration %ld ",i+1);
	for (j=0;j<dim;j++)
	  fprintf(file,"%e ",xdiff[i][j]);
	fprintf(file,"%e",diffcoeff[i]);
	fprintf(file,"\n");
      }
    }
    avpm=pm[0]*pm[0];
    loglikelihood= -log(pm[0]);
    for (i=1;i<dim;i++) {
      avpm += pm[i]*pm[i];
      loglikelihood -= log(pm[i]);
    }
    loglikelihood *= ((double)length);
    loglikelihood += -((double)length)*
      ((1.0+log(2.*M_PI))*dim)/2.0;
    avpm=sqrt(avpm/dim);
    fprintf(file,"#average forcast error= %e\n",avpm);
    fprintf(file,"#individual forecast errors: ");
    for (i=0;i<dim;i++)
      fprintf(file,"%e ",pm[i]);
    fprintf(file,"\n");
    if (arimaset)
      aic=2.0*(arpoles+mapoles)-2.0*loglikelihood;
    else
      aic=2.0*poles-2.0*loglikelihood;
    fprintf(file,"#Log-Likelihood= %e\t AIC= %e\n",loglikelihood,aic);
    for (i=0;i<size;i++) {
      id=aindex[0][i];
      is=aindex[1][i];
      if (id < dim)
	fprintf(file,"#x_%u(n-%u) ",id+1,is);
      else
	fprintf(file,"#e_%u(n-%u) ",id+1-dim,is);
      for (j=0;j<dim;j++)
	fprintf(file,"%e ",coeff[j][i]);
      fprintf(file,"\n");
    }
    if (!run_model || (verbosity&VER_USR1)) {
      for (i=poles;i<length;i++) {
	if (run_model)
	  fprintf(file,"#");
	for (j=0;j<dim;j++)
	  if (verbosity&VER_USR2)
	    fprintf(file,"%e %e ",series[j][i]+my_average[j],diff[j][i]);
	  else
	    fprintf(file,"%e ",diff[j][i]);
	fprintf(file,"\n");
      }
    }
    if (run_model && (ilength > 0)) {
      if (!arimaset)
	iterate_model(coeff,pm,diff,file);
      else
	iterate_arima_model(coeff,pm,diff,file);
    }
    fclose(file);
  }
  if (outfile != NULL)
    free(outfile);
  if (infile != NULL)
    free(infile);
  for (i=0;i<dim;i++) {
    free(coeff[i]);
    free(diff[i]);
    free(series[i]);
    if (arimaset)
      free(series[i+dim]);
  }
  free(coeff);
  free(diff);
  free(series);

  free(pm);

  return 0;
}
