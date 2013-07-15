/*=================================================================
 *
 * COMPLEXITYBS.C	.MEX file corresponding to COMPLEXITYBS.M
 *              returns the Lempel-Ziv compexity, and approximate 
 *              entropy for a n-bin encoding of the bit stream x
 *
 * The calling syntax is:
 *
 *		cmp = complexity(x)
 *
 *
 * This is a MEX-file for MATLAB.  
 *=================================================================*/
/* $Revision: 1.5 $ */
#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <math.h>

void complexity(double data[], const int length, double *cmp);
bool issubstring(unsigned short s[], const int ns,const int nq);
int compare(const void *arg1, const void *arg2);

void complexity(double	*data,
		int length,
		double *cmp)
     /*actually calculate the complexity*/
{
	unsigned short *s;
	double *sorted;
	float *td;
	int i,j,pc=0,nc=0;
	float mean=0,min=0,max=0;
	unsigned long int c=1;
	int ns,nq,k,bins;

	/*allocate memory for symbol sequence*/
	s = (unsigned short *) calloc(length,sizeof(unsigned short));	

	/* do the encoding */
	bins=1;
	for (i=0; i<length; i++)
	  {
	    *(s+i)= (int)floor(data[i])+1;
	    bins= (bins<s[i])? (s[i]):bins;
	  }

	/* now calculate the complexity */
	ns=1;
	nq=1;
	k=2;	
	while (k<length)
	{
		if (issubstring(s,ns,nq))
		{
		  /* Q is a substring of SQpi */
			nq++;
		}
		else
		{
		  /* Q is a new symbol sequence */
		  c++; 
		  ns=ns+nq;
		  nq=1;
		}
		k++;
	}

	/* free memory allocated for s */
	free(s);
	
	/*normalized complexity */
	*cmp=(c*log(length))/(length*log(bins));

	/*now, calculate the entropy*/
	/*NOT!!!!*/
}

int compare(const void *arg1, const void *arg2)
     /* compare two doubles and return arg1-arg2 */
{
  return( *(int *)arg1 - *(int *)arg2 );
}

bool issubstring(unsigned short s[], const int ns, const int nq)
     /* determine if Q is a substring of SQpi.*/
{
	int k,i;
	bool same;
	for (i=0; i<=(ns-nq); i++)
	{
		same=1;
		for (k=0; k<nq; k++)
		{
			if (s[i+k]!=s[ns+k])
				same=0;
		}
		if (same)
		{
		  return(1);
		}
	}
	return(0);
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     /* the MATLAB mex wrapper function */     
{ 
  double *x,*bins,*cmp;
  double thisCmp;
    int mrows,ncols,i; 
    int lengthx,bin,nbins;
    bool warnThem;

    
    /* Check for proper number of arguments */

    
    if (nrhs > 2) { 
	mexErrMsgTxt("Too many input arguments."); 
    } 
    if (nrhs == 0) {
        mexErrMsgTxt("Insufficient input arguments.");
    } 
    if (nrhs < 2) {
      /* set bin=2 */
      nbins=1;
      bins = (double *) calloc(1,sizeof(double));	
      *bins =2;
    }  
    if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /*Assign a pointer to the input matrix*/
    if (mxIsLogical(prhs[0])) {
      mexErrMsgTxt("Input argument may not be Logical");
    } else {
    x = mxGetPr(prhs[0]);
    }

    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */     
    mrows = mxGetM(prhs[0]); 
    ncols = mxGetN(prhs[0]);
    /* check that x is a vector */
    if ((mrows!=1) && (ncols!=1)) {
      mexWarnMsgTxt("First input should be a vector");
    }
    lengthx = mrows*ncols;
    
    /* Get the size of the "forbidden zone" --- the second input arg. */
    if (nrhs==2){
	mexWarnMsgTxt("Second input ignored");
    }

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); 

    /* Assign pointer to the ouput matrix */ 
    cmp = mxGetPr(plhs[0]);

    /* Do the actual computations in a subroutine */
    complexity(x,lengthx,&thisCmp);
    *cmp=thisCmp;
    

    return;
    
}


