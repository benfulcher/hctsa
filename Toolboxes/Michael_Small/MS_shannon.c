/*=================================================================
 *
 * SHANNON.C	.MEX file corresponding to SHANNON.M
 *              returns the approximate Shannon Entropy, for a n-bin
 *              encoding of the time series x at depth d
 *                 I.e. -\sum Plog(P)
 *              where P=(x_i,x_{i+1},....x_{i+d}) and the sum is
 *                              over all trajectories P
 *
 * The calling syntax is:
 *
 *		ent= entropy(x,n,d)
 *
 *
 * This is a MEX-file for MATLAB.
 *=================================================================*/
/* $Revision: 1.5 $ */
#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <math.h>

int compare(const void *arg1, const void *arg2)
     /* compare two doubles and return arg1-arg2 */
{
  return( *(int *)arg1 - *(int *)arg2 );
}

void entropy(double	*data,
		int length,
		int bin,
                int depth,
		double *ent)
     /*actually calculate the entropy*/
{
	unsigned short *s;
	double *sorted;
	float *td;
        unsigned long int *tally;
	int i,j,pc=0,nc=0;
	float total=0,prob=0;
	int k;
	unsigned long int os;

	/*allocate memory for symbol sequence*/
	s = (unsigned short *) calloc(length,sizeof(unsigned short));
	sorted = (double *) calloc(length,sizeof(double));
	td = (float *) calloc(bin-1,sizeof(float));

        /*determine threshold --- for n-bit encoding */
        /* copy the data */
        for (i=0; i<length; i++)
          *(sorted+i) = *(data+i);
        /*sort the copy*/
        qsort(sorted,length,sizeof(double),compare);
        /*extract the relevant percentiles*/
	for (i=1; i<bin; i++)
          *(td+i-1) = *(sorted+(int)(i*length/bin));

        /*dispose of the sorted data*/
        free(sorted);

	/* do the encoding */
	for (i=0; i<length; i++)
	  {
	    *(s+i)=0;
	    for (j=0; j<(bin-1); j++)
	      *(s+i) += (data[i]<*(td+j));
	  }

	/*allocate memory for the counters*/
	tally = (unsigned long int *) calloc(pow(bin,depth),sizeof(unsigned long int));
        for (i=0; i<pow(bin,depth); i++)
	  *(tally+i)=0;

	/* now calculate the entropy */
	length=length-depth+1;
	for (k=0; k<length; k++)
	{
	  os=0;
	  /*whats the encoded symbol sequence?*/
	  for (i=0; i<depth; i++)
	    os += pow(bin,i)*(*(s+i+k));
	  /* increase the appropriate counter */
	  (*(tally+os)) ++;
	}

	/*now calculate sum of P log P */
	total=0;
	for (i=0; i<pow(bin,depth); i++)
	  {
	    prob = (*(tally+i)) / (double)length; /* the probabilty of the i-th encoding - P */
	    if (prob>0)
	      total += prob*log(prob); /*the running sum of P*log(P) */
	  }

	/* free memory allocated for s */
	free(s);
	free(td);

	/* free memory allocated for tally */
	free(tally);

	/*now, calculate the entropy*/
	*ent=-total;

}

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )
     /* the MATLAB mex wrapper function */
{
    double *x,*bins,*deps,*ent;
    double thisEnt;
    int mrows,ncols,i,j;
    int lengthx,bin,nbins,dep,ndeps;
    bool warnThem;

    /* Check for proper number of arguments */

    if (nrhs > 3) {
	mexErrMsgTxt("Too many input arguments.");
    }
    if (nrhs == 0) {
        mexErrMsgTxt("Insufficient input arguments.");
    }
    if (nrhs < 2) {
      /* set bin=2 */
      nbins=1;
      bins = (double *) calloc(1,sizeof(double));
      *bins=2;
    }
    if (nrhs < 3) {
      /* set dep=3 */
      ndeps=1;
      deps = (double *) calloc(1,sizeof(double));
      *deps=3;
    }
    if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments.");
    }

    /*Assign a pointer to the input matrix*/
    x = mxGetPr(prhs[0]);

    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    /* check that x is a vector */
    if ((mrows!=1) && (ncols!=1)) {
      mexWarnMsgTxt("First input should be a vector");
    }
    lengthx = mrows*ncols;

    /* Get the size of the "forbidden zone" --- the second input arg. */
    if (nrhs>=2){
      bins=mxGetPr(prhs[1]);
      mrows = mxGetM(prhs[1]);
      ncols = mxGetN(prhs[1]);
      /* check that x is a vector */
      if ((mrows!=1) && (ncols!=1)) {
	mexWarnMsgTxt("Second input should be a vector");
      }
      nbins = mrows*ncols;
      /* check that all bins are integers */
      warnThem=0;
      for (i=0; i<nbins; i++) {
	if (floor(*(bins+i))!=*(bins+i)) {
	  warnThem=1;
	  *(bins+i)=(int)floor(*(bins+i));
	}
      }
      if (warnThem) {
	mexWarnMsgTxt("Second input being rounded to integer value(s).");
      }
    }

    /* */
    if (nrhs==3){
      deps=mxGetPr(prhs[2]);
      mrows = mxGetM(prhs[2]);
      ncols = mxGetN(prhs[2]);
      /* check that x is a vector */
      if ((mrows!=1) && (ncols!=1)) {
	mexWarnMsgTxt("Third input should be a vector");
      }
      ndeps = mrows*ncols;
      /* check that all bins are integers */
      warnThem=0;
      for (i=0; i<ndeps; i++) {
	if (floor(*(deps+i))!=*(deps+i)) {
	  warnThem=1;
	  *(deps+i)=(int)floor(*(deps+i));
	}
      }
      if (warnThem) {
	mexWarnMsgTxt("Third input being rounded to integer value(s).");
      }
    }


    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(ndeps, nbins, mxREAL);

    /* Assign pointer to the ouput matrix */
    ent = mxGetPr(plhs[0]);

    /* Do the actual computations in a subroutine */
    for (i=0; i<nbins; i++)
      for (j=0; j<ndeps; j++)
	{
	  bin=*(bins+i);
	  dep=*(deps+j);
	  entropy(x,lengthx,bin,dep,&thisEnt);
	  *(ent+i*ndeps+j)=thisEnt;
	}

    /*free allocated memory (if any) */
    if (nrhs<2)
      free(bins);
    if (nrhs<3)
      free(deps);

    return;

}
